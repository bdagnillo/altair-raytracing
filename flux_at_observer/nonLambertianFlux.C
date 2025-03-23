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

class BRDF {
private:
    double m_roughness; // Surface roughness parameter (0 = perfect mirror, 1 = fully diffuse)
    double m_specular;  // Specular reflection coefficient
    double m_diffuse;   // Diffuse reflection coefficient

public:
    BRDF(double roughness = 0.5, double specular = 0.3, double diffuse = 0.7) 
        : m_roughness(roughness), m_specular(specular), m_diffuse(diffuse) {
        // Ensure coefficients sum to 1
        double sum = m_specular + m_diffuse;
        m_specular /= sum;
        m_diffuse /= sum;
    }

    TVector3 SampleDirection(const TVector3& normal, const TVector3& incident) {
        // Decide between specular and diffuse reflection
        if (gRandom->Uniform() < m_specular) {
            return SampleSpecular(normal, incident);
        } else {
            return SampleDiffuse(normal);
        }
    }

private:
    TVector3 SampleSpecular(const TVector3& normal, const TVector3& incident) {
        // Perfect mirror reflection direction
        TVector3 reflect = incident - 2 * incident.Dot(normal) * normal;
        reflect.SetMag(1.0);

        // Add perturbation based on roughness
        double theta = gRandom->Gaus(0, m_roughness * M_PI/6); // Max 30° deviation at 1.0 roughness
        double phi = gRandom->Uniform(0, 2*M_PI);

        // Create perturbation vector perpendicular to reflection
        TVector3 perp1 = reflect.Orthogonal();
        TVector3 perp2 = reflect.Cross(perp1);
        
        // Apply perturbation
        TVector3 result = reflect + sin(theta)*(cos(phi)*perp1 + sin(phi)*perp2);
        result.SetMag(1.0);
        return result;
    }

    TVector3 SampleDiffuse(const TVector3& normal) {
        // Sample cosine-weighted hemisphere distribution
        double theta = acos(sqrt(gRandom->Uniform()));  // Cosine-weighted
        double phi = gRandom->Uniform(0, 2*M_PI);
        
        // Create local coordinate system
        TVector3 u = normal.Orthogonal();
        TVector3 v = normal.Cross(u);
        
        // Convert spherical to cartesian coordinates in local system
        double x = sin(theta) * cos(phi);
        double y = sin(theta) * sin(phi);
        double z = cos(theta);
        
        // Transform to world coordinates
        return (x * u + y * v + z * normal).Unit();
    }
};

// Add this to the global scope, before the traceRays function
BRDF gBRDF(0.3, 0.4, 0.6); // Default parameters for surface properties

void setupOpticsManager(AOpticsManager* manager) {
    manager->SetLimit(10000);
    TGeoBBox* box = new TGeoBBox("box", 200*cm, 200*cm, 200*cm);
    AOpticalComponent* world = new AOpticalComponent("world", box);
    manager->SetTopVolume(world);
    TGeoSphere* sphere = new TGeoSphere("sphereWithExitPort", 100.1*cm, 101*cm, 0., thetaMax);
    AMirror* mirror = new AMirror("mirror", sphere);
    ABorderSurfaceCondition* condition = new ABorderSurfaceCondition(world, mirror);
    condition->EnableLambertian(true);
    condition->SetGaussianRoughness(0.5); // Set roughness parameter
    world->AddNode(mirror, 1);
    manager->SetNsegments(100);
    manager->CloseGeometry();
}

bool isRayPassingThroughExitPort(ARay* ray, double exitPortZ) {
    Double_t lastPoint[3];
    ray->GetLastPoint(lastPoint);
    return lastPoint[2] < exitPortZ;
}

// Modify the traceRays function
int traceRays(AOpticsManager* manager, int n, double exitPortZ, Detector& detector, bool drawRays = false) {
    std::vector<TPolyLine3D*> rays;
    int hitCount = 0;
    
    // Add random number generator initialization
    gRandom->SetSeed();
    
    for (int i = 0; i < n; ++i) {
        double x = -60*cm;
        ARay* ray = new ARay(0, 400*nm, x, 0*cm, -80*cm, 0, 5, 0, 0);
        
        // Get initial direction for later use
        Double_t initDir[3];
        ray->GetDirection(initDir);
        TVector3 incidentDir(initDir[0], initDir[1], initDir[2]);
        
        manager->TraceNonSequential(*ray);
        
        // Get final position and normal at reflection point
        Double_t lastPoint[3];
        ray->GetLastPoint(lastPoint);
        TVector3 reflectionPoint(lastPoint[0], lastPoint[1], lastPoint[2]);
        
        // Calculate surface normal at reflection point (simplified - assumes spherical mirror)
        TVector3 normal = reflectionPoint.Unit();
        
        // Generate new direction using BRDF
        TVector3 newDir = gBRDF.SampleDirection(normal, incidentDir);
        
        // Create new ray from reflection point with new direction
        ARay* scatteredRay = new ARay(0, 400*nm, 
                                     reflectionPoint.X(), reflectionPoint.Y(), reflectionPoint.Z(),
                                     0, newDir.X(), newDir.Y(), newDir.Z());
        manager->TraceNonSequential(*scatteredRay);
        
        if (drawRays) {
            // Draw both initial and scattered rays
            TPolyLine3D* pol1 = ray->MakePolyLine3D();
            TPolyLine3D* pol2 = scatteredRay->MakePolyLine3D();
            
            if (isRayPassingThroughExitPort(scatteredRay, exitPortZ) && 
                detector.checkIntersection(scatteredRay)) {
                pol1->SetLineColor(kBlue);   // Initial ray that leads to detection
                pol2->SetLineColor(kGreen);  // Scattered ray that hits detector
                detector.hitCount++;
                hitCount++;
            } else if (isRayPassingThroughExitPort(scatteredRay, exitPortZ)) {
                pol1->SetLineColor(kGray);   // Initial ray
                pol2->SetLineColor(kYellow); // Scattered ray through exit port
            } else {
                pol1->SetLineColor(kGray);   // Initial ray
                pol2->SetLineColor(kRed);    // Scattered ray misses
            }
            pol1->SetLineWidth(1);
            pol2->SetLineWidth(1);
            pol1->Draw();
            pol2->Draw();
            rays.push_back(pol1);
            rays.push_back(pol2);
        } else if (isRayPassingThroughExitPort(scatteredRay, exitPortZ) && 
                   detector.checkIntersection(scatteredRay)) {
            detector.hitCount++;
            hitCount++;
        }

        delete ray;
        delete scatteredRay;
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
    
    std::cout << "\nFlux map data saved to 'fluxmap_data.csv'" << std::endl;
}

void visualizeDetector(double theta = 45.0, double phi = 0.0) {
    // Create canvas with square proportions
    TCanvas* c = new TCanvas("cvis", "Detector Visualization", 800, 800);
    c->cd();
    
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    setupOpticsManager(manager);
    
    // Create and position detector first
    Detector detector(20*cm, 20*cm);
    detector.setPosition(theta, phi, 100*cm);
    
    // Get the world volume and configure it
    TGeoVolume* world = gGeoManager->GetTopVolume();
    world->SetLineColor(kGray);
    world->SetTransparency(80);
    
    // Add detector to geometry
    detector.AddToGeometry(world);
    
    // Initialize geometry first
    gGeoManager->CloseGeometry();
    
    // Set up the 3D viewer
    gStyle->SetCanvasPreferGL(true);
    world->Draw("ogl");
    
    // Get viewer and set initial properties
    TGLViewer* glv = (TGLViewer*)gPad->GetViewer3D();
    if (glv) {
        glv->SetStyle(TGLRnrCtx::kWireFrame);
        glv->SetDrawOption(""); 
    }
    
    // Draw rays after geometry is initialized
    int n = 1000;  // Increased for better visualization
    double exitPortZ = -100*cm;
    int hitCount = traceRays(manager, n, exitPortZ, detector, true);
    
    // Final viewer updates
    if (glv) {
        glv->DrawGuides();
        glv->UpdateScene();
        glv->ResetCamerasAfterNextUpdate();
    }
    
    c->Update();
    c->Modified();
    gPad->Modified();
    gPad->Update();
    
    // Print detector information
    cout << "\nDetector Information:" << endl;
    cout << "Position (x,y,z): (" << detector.x/cm << ", " << detector.y/cm << ", " << detector.z/cm << ") cm" << endl;
    cout << "Normal vector (nx,ny,nz): (" << detector.nx << ", " << detector.ny << ", " << detector.nz << ")" << endl;
    cout << "Angular position: theta = " << theta << "°, phi = " << phi << "°" << endl;
    cout << "Number of rays hitting the detector: " << hitCount << endl;
}
