#include "TCanvas.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoSphere.h"
#include "TGeoTube.h"
#include "TH2.h"
#include "TThread.h"
#include "TVector3.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TGeoVolume.h"
#include "TGLViewer.h"
#include "TPad.h"

#include "AFocalSurface.h"
#include "AGeoUtil.h"
#include "AMirror.h"
#include "AObscuration.h"
#include "AOpticsManager.h"
#include "ABorderSurfaceCondition.h"
#include "ARayShooter.h"
#include <iostream>
#include <cmath>
#include <iomanip>

const double cm = AOpticsManager::cm();
const double nm = AOpticsManager::nm();
const double thetaMax = 170.;     // governs size of exit port

double calculateScatteringProbability(const TVector3& direction, const TVector3& normal) {
    double cosTheta = direction.Dot(normal);
    const double n = 2.0;  // Reduced from 5.0 to make scattering more visible
    double probability = pow(std::abs(cosTheta), n);
    return probability;
}

TVector3 generateScatteredDirection(const TVector3& incident, const TVector3& normal) {
    // Use a wider angle for more visible scattering
    const double maxAngle = 60.0 * TMath::Pi() / 180.0; // Increased from 45° to 60°
    
    // Create a local coordinate system based on the normal
    TVector3 w = normal.Unit();
    TVector3 u = TVector3(0, 1, 0).Cross(w).Unit();
    TVector3 v = w.Cross(u);
    
    while (true) {
        // Generate random angles with cosine weighting
        double r1 = (double)rand() / RAND_MAX;
        double r2 = (double)rand() / RAND_MAX;
        
        double theta = maxAngle * r1;  // Linear distribution in theta
        double phi = 2.0 * TMath::Pi() * r2;
        
        // Convert to cartesian coordinates
        double x = sin(theta) * cos(phi);
        double y = sin(theta) * sin(phi);
        double z = cos(theta);
        
        // Transform to world space
        TVector3 scattered = x*u + y*v + z*w;
        scattered = scattered.Unit();
        
        // Accept/reject based on our probability distribution
        double p = calculateScatteringProbability(scattered, normal);
        if ((double)rand() / RAND_MAX <= p) {
            return scattered;
        }
    }
}

struct Detector {
    double x, y, z;           // position
    double nx, ny, nz;        // normal vector
    double width, height;     // dimensions
    int hitCount;             // number of rays that hit the detector
    
    Detector(double w = 10*cm, double h = 10*cm) 
        : width(w), height(h), hitCount(0) {}
    
    /* Sets the detector position and orientation using spherical coordinates
     * @param theta: polar angle from negative z-axis (in degrees)
     *              0° = on axis with exit port
     *              90° = perpendicular to exit port axis
     * @param phi: azimuthal angle in x-y plane (in degrees)
     *            0° = +x axis
     *            90° = +y axis
     * @param radius: distance from exit port center
     */
    void setPosition(double theta, double phi, double radius) {
        // Convert spherical to cartesian coordinates
        double theta_rad = theta * M_PI / 180.0;
        double phi_rad = phi * M_PI / 180.0;
        
        // Calculate position relative to world origin
        x = radius * sin(theta_rad) * cos(phi_rad);
        y = radius * sin(theta_rad) * sin(phi_rad);
        z = -100*cm - radius * cos(theta_rad);  // Offset by exit port position
        
        // Calculate normal vector pointing to (0,0,-100)
        double dx = x - 0;
        double dy = y - 0;
        double dz = z - (-100*cm);
        double mag = sqrt(dx*dx + dy*dy + dz*dz);
        
        nx = -dy/mag;  // Negative because we want to point towards the exit port
        ny = dx/mag;
        nz = dz/mag;
    }
    
    bool checkIntersection(const ARay* ray) {
        Double_t lastPoint[3];
        Double_t direction[3];
        ray->GetLastPoint(lastPoint);
        ray->GetDirection(direction);
        
        // Calculate dot product of ray direction with detector normal
        double dot = direction[0]*nx + direction[1]*ny + direction[2]*nz;
        if (fabs(dot) < 1e-10) return false;  // Ray is parallel to detector plane
        
        // Calculate intersection with detector plane
        double dx = lastPoint[0] - x;  // Changed sign: point - center
        double dy = lastPoint[1] - y;
        double dz = lastPoint[2] - z;
        
        // Distance along ray to intersection point
        double t = -(dx*nx + dy*ny + dz*nz) / dot;  // Added negative sign
        
        // Calculate intersection point
        double ix = lastPoint[0] + direction[0] * t;
        double iy = lastPoint[1] + direction[1] * t;
        double iz = lastPoint[2] + direction[2] * t;
        
        // Vector from detector center to intersection point
        double rx = ix - x;
        double ry = iy - y;
        double rz = iz - z;
        
        // Compute distance in detector plane
        // Using cross product with normal to get perpendicular components
        double ux = ny*rz - nz*ry;  // First perpendicular direction
        double uy = nz*rx - nx*rz;
        double uz = nx*ry - ny*rx;
        
        double r2 = ux*ux + uy*uy + uz*uz;  // Square of distance in detector plane
        
        return r2 <= (width/2)*(width/2);
    }

    TGeoVolume* CreateGeometry() const {
        // Create a disk instead of a box
        TGeoTube* detDisk = new TGeoTube("detDisk", 
            0,           // inner radius
            width/2,     // outer radius (using width as diameter)
            1*cm         // half-length in z (thickness)
        );
        TGeoVolume* detVol = new TGeoVolume("detector", detDisk);
        detVol->SetLineColor(kBlue);
        detVol->SetLineWidth(2);
        detVol->SetFillColor(kCyan);
        detVol->SetTransparency(20);
        return detVol;
    }
    
    void AddToGeometry(TGeoVolume* top) const {
        TGeoVolume* detVol = CreateGeometry();
        
        // Calculate rotation matrix to point detector normal to exit port
        // We need to rotate from the default orientation (normal along z-axis)
        // to point towards (0,0,-100)
        
        // First convert normal vector to spherical coordinates
        double theta_rot = acos(-nz) * 180/M_PI;  // Note the negative nz
        double phi_rot = atan2(-ny, -nx) * 180/M_PI;  // Note the negatives
        
        // Create rotation that first rotates around Y axis (theta)
        // then around Z axis (phi)
        TGeoRotation* rot = new TGeoRotation();
        rot->SetAngles(phi_rot, theta_rot, 0);
        
        TGeoTranslation* trans = new TGeoTranslation(x, y, z);
        TGeoCombiTrans* transform = new TGeoCombiTrans(*trans, *rot);
        
        top->AddNode(detVol, 1, transform);
    }
};

class NonLambertianSurface : public ABorderSurfaceCondition {
public:
    NonLambertianSurface(AOpticalComponent* comp1, AOpticalComponent* comp2) 
        : ABorderSurfaceCondition(comp1, comp2) {}
    
    virtual bool Reflection(ARay& ray, const TVector3& normal) {
        // Debug output
        std::cout << "NonLambertianSurface::Reflection called" << std::endl;
        
        // Get incident direction
        Double_t dir[3];
        ray.GetDirection(dir);
        TVector3 incident(dir[0], dir[1], dir[2]);
        
        // Print incident and normal directions
        std::cout << "Incident: (" << dir[0] << ", " << dir[1] << ", " << dir[2] << ")" << std::endl;
        std::cout << "Normal: (" << normal.X() << ", " << normal.Y() << ", " << normal.Z() << ")" << std::endl;
        
        // Generate scattered direction
        TVector3 scattered = generateScatteredDirection(incident, normal);
        
        // Ensure the scattered direction is in the correct hemisphere
        if (scattered.Dot(normal) < 0) {
            scattered = -scattered;
        }
        
        // Print scattered direction
        std::cout << "Scattered: (" << scattered.X() << ", " << scattered.Y() << ", " << scattered.Z() << ")" << std::endl;
        
        // Set the new direction
        ray.SetDirection(scattered.X(), scattered.Y(), scattered.Z());
        return true;
    }
};

void setupOpticsManager(AOpticsManager* manager) {
    manager->SetLimit(10000);
    
    // Create world volume with transparent vacuum
    TGeoBBox* box = new TGeoBBox("box", 200*cm, 200*cm, 200*cm);
    AOpticalComponent* world = new AOpticalComponent("world", box);
    world->SetMedium(world->GetTransparentVacuumMedium());
    manager->SetTopVolume(world);
    
    // Create mirror with opaque vacuum
    TGeoSphere* sphere = new TGeoSphere("sphereWithExitPort", 100.1*cm, 101*cm, 0., thetaMax);
    AMirror* mirror = new AMirror("mirror", sphere);
    mirror->SetMedium(world->GetOpaqueVacuumMedium());
    
    // Create and set the non-Lambertian surface condition
    NonLambertianSurface* condition = new NonLambertianSurface(world, mirror);
    condition->EnableLambertian(false);
    condition->SetGaussianRoughness(0.5); // Set roughness parameter
    world->AddBorderSurfaceCondition(condition);
    
    // Add mirror to world
    world->AddNode(mirror, 1);
    
    // Set geometry parameters
    manager->SetNsegments(100);
    
    // Explicitly create voxels BEFORE closing geometry
    // AOpticalComponent inherits from TGeoVolume, so it already is a volume
    world->Voxelize("");  // Force voxelization with empty option
    
    // Now close geometry
    manager->CloseGeometry();
    
    // Verify voxelization was successful
    TGeoVoxelFinder* voxels = world->GetVoxels();
    if (!voxels) {
        std::cerr << "WARNING: Voxelization failed!" << std::endl;
    }
}

bool isRayPassingThroughExitPort(ARay* ray, double exitPortZ) {
    Double_t lastPoint[3];
    ray->GetLastPoint(lastPoint);
    return lastPoint[2] < exitPortZ;
}

// Fix traceRays function to update detector's hit count
int traceRays(AOpticsManager* manager, int n, double exitPortZ, Detector& detector, bool drawRays = false) {
    std::vector<TPolyLine3D*> rays;  // Store rays for visualization
    int hitCount = 0;  // Count of rays hitting the detector
    
    for (int i = 0; i < n; ++i) {
        double x = -60*cm;
        ARay* ray = new ARay(0, 400*nm, x, 0*cm, -80*cm, 0, 5, 0, 0);
        manager->TraceNonSequential(*ray);
        
        if (drawRays) {
            TPolyLine3D* pol = ray->MakePolyLine3D();
            if (isRayPassingThroughExitPort(ray, exitPortZ) && 
                detector.checkIntersection(ray)) {
                pol->SetLineColor(kGreen);  // Hits detector
                detector.hitCount++;  // Update detector's counter
                hitCount++;
            } else if (isRayPassingThroughExitPort(ray, exitPortZ)) {
                pol->SetLineColor(kYellow);  // Passes through exit port but misses detector
            } else {
                pol->SetLineColor(kRed);    // Doesn't pass through exit port
            }
            pol->SetLineWidth(2);
            pol->Draw();
            rays.push_back(pol);  // Store for later cleanup if needed
        } else if (isRayPassingThroughExitPort(ray, exitPortZ) && 
                   detector.checkIntersection(ray)) {
            detector.hitCount++;  // Update detector's counter
            hitCount++;
        }

        delete ray;
    }
    return hitCount;
}

// Fix fraction calculation in sweepDetector
void sweepDetector() {
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    setupOpticsManager(manager);

    int n = 100000; // reduced number of rays per point for faster sweep
    double exitPortZ = -100*cm;
    
    // Create 2D histogram for the angular sweep
    const int nThetaBins = 45;  // 2-degree steps in theta
    const int nPhiBins = 20;    // 5-degree steps in phi
    TH2D* fluxMap = new TH2D("fluxMap", "Detector Flux Map;#theta (deg);#phi (deg)", 
                            nThetaBins, 0, 90, nPhiBins, 0, 360);
    
    // Create detector
    Detector detector;
    
    // Perform the sweep
    std::cout << "\nStarting detector sweep..." << std::endl;
    std::cout << "Format: theta(°), phi(°): hits/total = fraction" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    
    for(int i = 0; i < nThetaBins; i++) {
        double theta = (i + 0.5) * 90.0/nThetaBins;
        for(int j = 0; j < nPhiBins; j++) {
            double phi = (j + 0.5) * 360.0/nPhiBins;
            
            detector.hitCount = 0; // reset counter
            detector.setPosition(theta, phi, 100*cm);
            
            traceRays(manager, n, exitPortZ, detector);
            
            // Fix fraction calculation
            double fraction = double(detector.hitCount)/double(n);
            fluxMap->SetBinContent(i+1, j+1, fraction);
            
            // Detailed progress update
            std::cout << std::fixed << std::setprecision(1)
                     << theta << "°, " << phi << "°: "
                     << detector.hitCount << "/" << n << " = "
                     << std::setprecision(8) << fraction
                     << std::endl;
        }
        // Add separator between theta values
        std::cout << "----------------------------------------" << std::endl;
    }
    
    // Create canvas and draw the histogram
    TCanvas* c = new TCanvas("c", "Detector Flux Map", 1000, 800);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);
    
    fluxMap->GetZaxis()->SetTitle("Fraction of rays detected");
    fluxMap->Draw("COLZ");
    
    // Add color scale
    gPad->Update();
    TPaletteAxis* palette = (TPaletteAxis*)fluxMap->GetListOfFunctions()->FindObject("palette");
    if(palette) {
        palette->SetX1NDC(0.9);
        palette->SetX2NDC(0.92);
    }
    
    c->Update();
    
    // Save flux map data to CSV
    std::ofstream csvFile("fluxmap_data.csv");
    csvFile << "theta,phi,fraction\n";  // Header
    
    for(int i = 0; i < nThetaBins; i++) {
        double theta = (i + 0.5) * 90.0/nThetaBins;
        for(int j = 0; j < nPhiBins; j++) {
            double phi = (j + 0.5) * 360.0/nPhiBins;
            double fraction = fluxMap->GetBinContent(i+1, j+1);
            csvFile << std::fixed << std::setprecision(6)
                   << theta << "," << phi << "," << fraction << "\n";
        }
    }
    csvFile.close();
    
    // std::cout << "\nFlux map data saved to 'fluxmap_data.csv'" << std::endl;
}

void visualizeDetector(double theta = 45.0, double phi = 0.0) {
    // Create canvas with square proportions
    TCanvas* c = new TCanvas("cvis", "Detector Visualization", 800, 800);
    c->cd();
    
    // Delete any existing geometry manager first
    if (gGeoManager) {
        delete gGeoManager;
        gGeoManager = nullptr;
    }
    
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    setupOpticsManager(manager);
    
    // Create and position detector
    Detector detector(20*cm, 20*cm);
    detector.setPosition(theta, phi, 100*cm);

    // --------- RAYTRACING WITHOUT VISUALIZATION ----------
    // Trace rays without trying to visualize them
    int n = 1000;  // More rays for better statistics
    double exitPortZ = -100*cm;
    int hitCount = traceRays(manager, n, exitPortZ, detector, false);
    
    // --------- 2D VISUALIZATION INSTEAD OF 3D ----------
    // Draw a 2D representation of the system
    
    // Create a top view
    TPad* pad1 = new TPad("pad1", "Top View", 0, 0.5, 1, 1);
    pad1->Draw();
    pad1->cd();
    
    // Draw mirror outline (top view)
    TEllipse* mirror = new TEllipse(0, 0, 100*cm, 100*cm);
    mirror->SetFillStyle(0);
    mirror->SetLineColor(kBlue);
    mirror->SetLineWidth(2);
    mirror->Draw();
    
    // Draw exit port
    double exitPortRadius = 100*cm * sin(thetaMax * TMath::Pi() / 180.0);
    TEllipse* exitPort = new TEllipse(0, -100*cm, exitPortRadius, exitPortRadius);
    exitPort->SetFillStyle(0);
    exitPort->SetLineColor(kRed);
    exitPort->SetLineWidth(2);
    exitPort->Draw();
    
    // Draw detector position
    TEllipse* detCircle = new TEllipse(detector.x, detector.z + 100*cm, detector.width/2, detector.width/2);
    detCircle->SetFillStyle(0);
    detCircle->SetLineColor(kGreen);
    detCircle->SetLineWidth(2);
    detCircle->Draw();
    
    // Draw detector normal
    TArrow* normal = new TArrow(detector.x, detector.z + 100*cm, 
                        detector.x + detector.nx * 20*cm, 
                        (detector.z + 100*cm) + detector.nz * 20*cm, 
                        0.01, "|>");
    normal->SetLineColor(kGreen);
    normal->SetLineWidth(2);
    normal->Draw();
    
    // Add axes
    TLine* xAxis = new TLine(-150*cm, 0, 150*cm, 0);
    TLine* zAxis = new TLine(0, -200*cm, 0, 100*cm);
    xAxis->SetLineStyle(2);
    zAxis->SetLineStyle(2);
    xAxis->Draw();
    zAxis->Draw();
    
    // Add labels
    TLatex* xLabel = new TLatex(150*cm, -10*cm, "x");
    TLatex* zLabel = new TLatex(5*cm, -200*cm, "z");
    xLabel->Draw();
    zLabel->Draw();
    TLatex* topView = new TLatex(-140*cm, 80*cm, "Top View (x-z plane)");
    topView->SetTextSize(0.05);
    topView->Draw();
    
    // Create a side view
    c->cd();
    TPad* pad2 = new TPad("pad2", "Side View", 0, 0, 1, 0.5);
    pad2->Draw();
    pad2->cd();
    
    // Draw mirror outline (side view)
    TEllipse* mirrorSide = new TEllipse(0, 0, 100*cm, 100*cm);
    mirrorSide->SetFillStyle(0);
    mirrorSide->SetLineColor(kBlue);
    mirrorSide->SetLineWidth(2);
    mirrorSide->Draw();
    
    // Draw exit port
    TEllipse* exitPortSide = new TEllipse(-100*cm, 0, exitPortRadius, exitPortRadius);
    exitPortSide->SetFillStyle(0);
    exitPortSide->SetLineColor(kRed);
    exitPortSide->SetLineWidth(2);
    exitPortSide->Draw();
    
    // Draw detector position
    TEllipse* detCircleSide = new TEllipse(detector.z + 100*cm, detector.y, detector.width/2, detector.width/2);
    detCircleSide->SetFillStyle(0);
    detCircleSide->SetLineColor(kGreen);
    detCircleSide->SetLineWidth(2);
    detCircleSide->Draw();
    
    // Draw detector normal
    TArrow* normalSide = new TArrow(detector.z + 100*cm, detector.y, 
                         (detector.z + 100*cm) + detector.nz * 20*cm, 
                         detector.y + detector.ny * 20*cm, 
                         0.01, "|>");
    normalSide->SetLineColor(kGreen);
    normalSide->SetLineWidth(2);
    normalSide->Draw();
    
    // Add axes
    TLine* zAxisSide = new TLine(-200*cm, 0, 100*cm, 0);
    TLine* yAxisSide = new TLine(0, -150*cm, 0, 150*cm);
    zAxisSide->SetLineStyle(2);
    yAxisSide->SetLineStyle(2);
    zAxisSide->Draw();
    yAxisSide->Draw();
    
    // Add labels
    TLatex* zLabelSide = new TLatex(-200*cm, 10*cm, "z");
    TLatex* yLabelSide = new TLatex(5*cm, 150*cm, "y");
    zLabelSide->Draw();
    yLabelSide->Draw();
    TLatex* sideView = new TLatex(80*cm, -120*cm, "Side View (z-y plane)");
    sideView->SetTextSize(0.05);
    sideView->Draw();
    
    // Add detector info textbox
    TPaveText* info = new TPaveText(20*cm, -80*cm, 100*cm, -20*cm);
    info->SetFillColor(kWhite);
    info->SetTextAlign(12);
    info->AddText(Form("Detector at #theta = %.1f#circ, #phi = %.1f#circ", theta, phi));
    info->AddText(Form("Position: (%.1f, %.1f, %.1f) cm", detector.x/cm, detector.y/cm, detector.z/cm));
    info->AddText(Form("Rays detected: %d/%d (%.1f%%)", hitCount, n, 100.0*hitCount/n));
    info->Draw();
    
    c->Update();
    
    // Print detector information to console
    cout << "\nDetector Information:" << endl;
    cout << "Position (x,y,z): (" << detector.x/cm << ", " << detector.y/cm << ", " << detector.z/cm << ") cm" << endl;
    cout << "Normal vector (nx,ny,nz): (" << detector.nx << ", " << detector.ny << ", " << detector.nz << ")" << endl;
    cout << "Angular position: theta = " << theta << "°, phi = " << phi << "°" << endl;
    cout << "Number of rays hitting the detector: " << hitCount << " of " << n << " (" << 100.0*hitCount/n << "%)" << endl;
}

void visualizeDetectorSimple(double theta = 45.0, double phi = 0.0) {
    TCanvas* c = new TCanvas("cvis", "Detector Visualization", 800, 800);
    c->Divide(1, 2);
    
    c->cd(1);
    
    // Create manager and detector
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    setupOpticsManager(manager);
    
    Detector detector(20*cm, 20*cm);
    detector.setPosition(theta, phi, 100*cm);
    
    // Draw geometry in top view
    TView* topView = TView::CreateView(1);
    topView->SetRange(-150, -150, -200, 150, 150, 0);
    
    // Fix SetView call by adding the required 4th argument
    Int_t irep;
    topView->SetView(0, 0, 0, irep);
    
    // Draw mirror outline
    TEllipse* mirrorOutline = new TEllipse(0, 0, 100*cm, 100*cm);
    mirrorOutline->SetLineColor(kBlue);
    mirrorOutline->SetLineWidth(2);
    mirrorOutline->SetFillStyle(0);
    mirrorOutline->Draw();
    
    // Draw detector position as a marker
    TMarker* detectorMarker = new TMarker(detector.x, detector.y, 20);
    detectorMarker->SetMarkerColor(kRed);
    detectorMarker->SetMarkerStyle(20);
    detectorMarker->SetMarkerSize(1.5);
    detectorMarker->Draw();
    
    // Draw axes
    TLine* xaxis = new TLine(-150, 0, 150, 0);
    TLine* yaxis = new TLine(0, -150, 0, 150);
    xaxis->SetLineColor(kBlack);
    yaxis->SetLineColor(kBlack);
    xaxis->Draw();
    yaxis->Draw();
    
    c->cd(2);
    
    // Statistics panel
    TPaveText* textBox = new TPaveText(0.05, 0.05, 0.95, 0.95);
    textBox->SetFillColor(kWhite);
    textBox->SetTextAlign(12);  // Align left
    
    int n = 1000;
    double exitPortZ = -100*cm;
    int hitCount = traceRays(manager, n, exitPortZ, detector, false);
    
    textBox->AddText(Form("Detector angle: theta = %.1f°, phi = %.1f°", theta, phi));
    textBox->AddText(Form("Position (x,y,z): (%.1f, %.1f, %.1f) cm", detector.x/cm, detector.y/cm, detector.z/cm));
    textBox->AddText(Form("Normal vector: (%.3f, %.3f, %.3f)", detector.nx, detector.ny, detector.nz));
    textBox->AddText(Form("Rays traced: %d", n));
    textBox->AddText(Form("Rays detected: %d (%.2f%%)", hitCount, 100.0*hitCount/n));
    
    textBox->Draw();
    c->Update();
}

void visualizeDetectorText(double theta = 45.0, double phi = 0.0) {
    // Create and setup manager
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    setupOpticsManager(manager);
    
    // Create detector
    Detector detector(20*cm, 20*cm);
    detector.setPosition(theta, phi, 100*cm);
    
    // Trace rays
    int n = 10000;  // More rays for better statistics
    double exitPortZ = -100*cm;
    int hitCount = traceRays(manager, n, exitPortZ, detector, false);
    
    // Print detailed information
    std::cout << "\n===================================" << std::endl;
    std::cout << "DETECTOR INFORMATION" << std::endl;
    std::cout << "===================================" << std::endl;
    std::cout << "Angular position: theta = " << theta << "°, phi = " << phi << "°" << std::endl;
    std::cout << "Position (x,y,z): (" 
              << detector.x/cm << ", " 
              << detector.y/cm << ", " 
              << detector.z/cm << ") cm" << std::endl;
    std::cout << "Normal vector: (" 
              << detector.nx << ", " 
              << detector.ny << ", " 
              << detector.nz << ")" << std::endl;
    std::cout << "Rays traced: " << n << std::endl;
    std::cout << "Rays detected: " << hitCount << " (" 
              << 100.0*hitCount/n << "%)" << std::endl;
    std::cout << "===================================" << std::endl;
    
    // Create text-based ASCII visualization
    std::cout << "\nTop View (X-Z plane, Y=0):" << std::endl;
    std::cout << "    ^Z" << std::endl;
    std::cout << "    |" << std::endl;
    std::cout << "    |   ,-------," << std::endl;
    std::cout << "    |  /         \\" << std::endl;
    std::cout << "    | |     *     | Mirror" << std::endl;
    std::cout << "    |  \\         /" << std::endl;
    std::cout << "    |   '-------'" << std::endl;
    std::cout << "    |" << std::endl;
    
    // Calculate detector position in character space
    int detX = (int)(detector.x/cm / 20) + 10;  // Scale and offset
    int detZ = -(int)((detector.z+100*cm)/cm / 10) + 15;  // Scale, flip and offset
    
    // Draw detector position
    for (int z = 15; z >= 0; z--) {
        std::cout << "    |";
        for (int x = 0; x < 20; x++) {
            if (z == detZ && x == detX) {
                std::cout << "D";  // Detector
            } else {
                std::cout << " ";
            }
        }
        if (z == 15) std::cout << "  Exit Port at Z=-100cm";
        if (z == detZ) std::cout << "  <- Detector";
        std::cout << std::endl;
    }
    std::cout << "    +---------------------> X" << std::endl;
}

// Add this function to help debug geometry issues
void debugGeometry() {
    // Delete any existing geometry manager first
    if (gGeoManager) {
        delete gGeoManager;
        gGeoManager = nullptr;
    }
    
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    
    std::cout << "Setting up optics manager..." << std::endl;
    manager->SetLimit(10000);
    
    // Create world volume with transparent vacuum
    TGeoBBox* box = new TGeoBBox("box", 200*cm, 200*cm, 200*cm);
    AOpticalComponent* world = new AOpticalComponent("world", box);
    world->SetMedium(world->GetTransparentVacuumMedium());
    manager->SetTopVolume(world);
    
    std::cout << "Created world volume" << std::endl;
    
    // Create mirror with opaque vacuum
    TGeoSphere* sphere = new TGeoSphere("sphereWithExitPort", 100.1*cm, 101*cm, 0., thetaMax);
    AMirror* mirror = new AMirror("mirror", sphere);
    mirror->SetMedium(world->GetOpaqueVacuumMedium());
    
    std::cout << "Created mirror" << std::endl;
    
    // Create and set the non-Lambertian surface condition
    NonLambertianSurface* condition = new NonLambertianSurface(world, mirror);
    condition->EnableLambertian(false);
    world->AddBorderSurfaceCondition(condition);
    
    std::cout << "Added border condition" << std::endl;
    
    // Add mirror to world
    world->AddNode(mirror, 1);
    
    std::cout << "Added mirror to world" << std::endl;
    
    // Set geometry parameters and close
    manager->SetNsegments(100);
    
    std::cout << "Closing geometry..." << std::endl;
    manager->CloseGeometry();
    std::cout << "Geometry closed" << std::endl;
    
    // Check voxelization
    if (gGeoManager) {
        std::cout << "Checking voxelization for world volume..." << std::endl;
        TGeoVolume* worldVol = gGeoManager->GetTopVolume();
        TGeoVoxelFinder* voxels = worldVol->GetVoxels();
        if (voxels) {
            std::cout << "World volume has voxels" << std::endl;
        } else {
            std::cout << "World volume does NOT have voxels" << std::endl;
            std::cout << "Attempting to voxelize..." << std::endl;
            // Fix Voxelize call by providing the required option parameter
            worldVol->Voxelize("");
        }
    } else {
        std::cout << "gGeoManager is null" << std::endl;
    }
}

void visualizeDetectorSafe(double theta = 45.0, double phi = 0.0) {
    // Create and setup manager
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    setupOpticsManager(manager);
    
    // Create detector
    Detector detector(20*cm, 20*cm);
    detector.setPosition(theta, phi, 100*cm);
    
    // Trace rays
    int n = 10000;
    double exitPortZ = -100*cm;
    int hitCount = traceRays(manager, n, exitPortZ, detector, false);
    
    // Create canvas for 2D visualization
    TCanvas* c = new TCanvas("cvis", "Detector Safe Visualization", 900, 600);
    c->Divide(2, 1);
    
    // Top view (X-Z plane)
    c->cd(1);
    
    // Set coordinate range
    double rangeX = 200*cm;
    double rangeZ = 200*cm;
    
    TH2F* frame1 = new TH2F("frame1", "Top View (X-Z);X (cm);Z (cm)", 
                           100, -rangeX, rangeX, 100, -rangeZ, rangeZ);
    frame1->SetStats(0);
    frame1->Draw();
    
    // Draw mirror
    TEllipse* mirror = new TEllipse(0, 0, 100*cm, 100*cm);
    mirror->SetFillStyle(0);
    mirror->SetLineColor(kBlue);
    mirror->SetLineWidth(2);
    mirror->Draw();
    
    // Draw exit port
    double exitPortRadius = 100*cm * sin(thetaMax * TMath::Pi() / 180.0);
    TEllipse* exitPort = new TEllipse(0, -100*cm, exitPortRadius, exitPortRadius);
    exitPort->SetFillStyle(0);
    exitPort->SetLineColor(kRed);
    exitPort->SetLineWidth(2);
    exitPort->Draw();
    
    // Draw detector position
    TEllipse* detCircle = new TEllipse(detector.x, detector.z + 100*cm, detector.width/2, detector.height/2);
    detCircle->SetFillStyle(0);
    detCircle->SetLineColor(kGreen);
    detCircle->SetLineWidth(2);
    detCircle->Draw();
    
    // Side view (Z-Y plane)
    c->cd(2);
    
    TH2F* frame2 = new TH2F("frame2", "Side View (Z-Y);Z (cm);Y (cm)", 
                           100, -rangeZ, rangeZ, 100, -rangeX, rangeX);
    frame2->SetStats(0);
    frame2->Draw();
    
    // Draw mirror (side view)
    TEllipse* mirrorSide = new TEllipse(0, 0, 100*cm, 100*cm);
    mirrorSide->SetFillStyle(0);
    mirrorSide->SetLineColor(kBlue);
    mirrorSide->SetLineWidth(2);
    mirrorSide->Draw();
    
    // Draw exit port (side view)
    TEllipse* exitPortSide = new TEllipse(-100*cm, 0, exitPortRadius, exitPortRadius);
    exitPortSide->SetFillStyle(0);
    exitPortSide->SetLineColor(kRed);
    exitPortSide->SetLineWidth(2);
    exitPortSide->Draw();
    
    // Draw detector position (side view)
    TEllipse* detCircleSide = new TEllipse(detector.z + 100*cm, detector.y, detector.width/2, detector.height/2);
    detCircleSide->SetFillStyle(0);
    detCircleSide->SetLineColor(kGreen);
    detCircleSide->SetLineWidth(2);
    detCircleSide->Draw();
    
    // Add text with detector information
    TPaveText* info = new TPaveText(0.1, 0.1, 0.9, 0.3, "NDC");
    info->SetFillColor(kWhite);
    info->SetTextAlign(12);
    info->AddText(Form("Detector at #theta = %.1f#circ, #phi = %.1f#circ", theta, phi));
    info->AddText(Form("Position: (%.1f, %.1f, %.1f) cm", detector.x/cm, detector.y/cm, detector.z/cm));
    info->AddText(Form("Rays detected: %d/%d (%.1f%%)", hitCount, n, 100.0*hitCount/n));
    info->Draw();
    
    c->Update();
    
    // Print detector information to console as well
    std::cout << "\nDetector Information:" << std::endl;
    std::cout << "Position (x,y,z): (" << detector.x/cm << ", " << detector.y/cm << ", " << detector.z/cm << ") cm" << std::endl;
    std::cout << "Normal vector (nx,ny,nz): (" << detector.nx << ", " << detector.ny << ", " << detector.nz << ")" << std::endl;
    std::cout << "Angular position: theta = " << theta << "°, phi = " << phi << "°" << std::endl;
    std::cout << "Number of rays hitting the detector: " << hitCount << " of " << n << " (" << 100.0*hitCount/n << "%)" << std::endl;
}
