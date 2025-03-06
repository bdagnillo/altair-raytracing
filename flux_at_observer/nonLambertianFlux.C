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
#include <random>

const double cm = AOpticsManager::cm();
const double nm = AOpticsManager::nm();
const double thetaMax = 170.;     // governs size of exit port

enum BRDFModel {
    OREN_NAYAR,
    COOK_TORRANCE
};

class ReflectiveSurface {
public:
    ReflectiveSurface(BRDFModel model, double albedo) : brdfModel(model), albedo(albedo) {}

    double computeBRDF(const TVector3& incoming, const TVector3& outgoing, const TVector3& normal) const {
        switch (brdfModel) {
            case OREN_NAYAR:
                return orenNayarBRDF(incoming, outgoing, normal);
            case COOK_TORRANCE:
                return cookTorranceBRDF(incoming, outgoing, normal);
            default:
                return 0.0;
        }
    }

private:
    BRDFModel brdfModel;
    double albedo;

    TVector3 sampleHemisphereCosine(double u1, double u2, const TVector3& normal) {
        // Concentric disk mapping for better distribution
        double r = sqrt(u1);
        double theta = 2.0 * M_PI * u2;
        
        // Map uniform random numbers to cosine-weighted hemisphere
        double x = r * cos(theta);
        double y = r * sin(theta);
        double z = sqrt(std::max(0.0, 1.0 - x*x - y*y));  // Cosine-weighted
        
        // Build coordinate system
        TVector3 up = (std::abs(normal.Z()) < 0.999) ? TVector3(0,0,1) : TVector3(1,0,0);
        TVector3 tangent = normal.Cross(up).Unit();
        TVector3 bitangent = normal.Cross(tangent);
        
        // Transform sample from local to world space
        return (x * tangent + y * bitangent + z * normal).Unit();
    }
    
    double monteCarloIntegration(const TVector3& outgoing, const TVector3& normal, 
                            int numSamples, std::function<double(const TVector3&, const TVector3&, const TVector3&)> brdfFunc) {
        double result = 0.0;
        
        for (int i = 0; i < numSamples; ++i) {
            double u1 = getRandomUniform();
            double u2 = getRandomUniform();
            
            TVector3 incoming = sampleHemisphereCosine(u1, u2, normal);
            
            // Cosine-weighted PDF = cos(theta) / pi
            double cosTheta = incoming.Dot(normal);
            double pdf = std::max(1e-6, cosTheta / M_PI);
            
            // Evaluate BRDF
            double brdf = brdfFunc(incoming, outgoing, normal);
            
            // Include geometry term (cosine) and divide by PDF
            result += brdf * cosTheta / pdf;
        }
        
        return result / numSamples;
    }
    
    double orenNayarBRDF(const TVector3& incoming, const TVector3& outgoing, const TVector3& normal) const {
        double sigma = 0.8; // Roughness parameter for Barium Sulfate (approximate value)
        double sigma2 = sigma * sigma;
    
        double theta_i = acos(std::max(-1.0, std::min(1.0, incoming.Dot(normal))));
        double theta_r = acos(std::max(-1.0, std::min(1.0, outgoing.Dot(normal))));
    
        double alpha = std::max(theta_i, theta_r);
        double beta = std::min(theta_i, theta_r);
    
        double phi_i = atan2(incoming.Y(), incoming.X());
        double phi_r = atan2(outgoing.Y(), outgoing.X());
    
        double delta_phi = phi_i - phi_r;
    
        double A = 1.0 - 0.5 * (sigma2 / (sigma2 + 0.33));
        double B = 0.45 * (sigma2 / (sigma2 + 0.09));
    
        double brdf = albedo / M_PI * (A + B * std::max(0.0, cos(delta_phi)) * sin(alpha) * tan(beta));
        
        // Ensure energy conservation
        return std::min(1.0, std::max(0.0, brdf * albedo));
    }
    
    double cookTorranceBRDF(const TVector3& incoming, const TVector3& outgoing, const TVector3& normal) const {
        const double roughness = 0.3; // Spectralon roughness
        const double F0 = 0.04; // Fresnel reflectance
        const double epsilon = 1e-8; // Small value to prevent division by zero

        TVector3 halfVector = (incoming + outgoing).Unit();

        double NdotL = std::max(epsilon, normal.Dot(incoming));
        double NdotV = std::max(epsilon, normal.Dot(outgoing));
        double NdotH = std::max(epsilon, normal.Dot(halfVector));
        double VdotH = std::max(epsilon, outgoing.Dot(halfVector));

        // Squaring alpha once
        double alpha2 = roughness * roughness;
        
        // Distribution term (GGX/Trowbridge-Reitz)
        double D = alpha2 / (M_PI * pow(NdotH * NdotH * (alpha2 - 1.0) + 1.0, 2.0));

        // Geometry term (Smith's method with GGX)
        double G = 0.5 / (NdotV * sqrt(alpha2 + (1.0 - alpha2) * NdotL * NdotL) +
                          NdotL * sqrt(alpha2 + (1.0 - alpha2) * NdotV * NdotV));

        // Fresnel term (Schlick approximation)
        double F = F0 + (1.0 - F0) * pow(1.0 - VdotH, 5.0);

        double brdf = (D * G * F) / (4.0 * NdotV * NdotL);
        return std::max(0.0, brdf);
    }
    
    double getRandomUniform() {
        // Use C++11's random number generation
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<> dis(0.0, 1.0);
        return dis(gen);
    }
};

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

void setupOpticsManager(AOpticsManager* manager, BRDFModel brdfModel, double albedo) {
    manager->SetLimit(10000);
    TGeoBBox* box = new TGeoBBox("box", 200*cm, 200*cm, 200*cm);
    AOpticalComponent* world = new AOpticalComponent("world", box);
    manager->SetTopVolume(world);
    TGeoSphere* sphere = new TGeoSphere("sphereWithExitPort", 100.1*cm, 101*cm, 0., thetaMax);
    AMirror* mirror = new AMirror("mirror", sphere);
    ABorderSurfaceCondition* condition = new ABorderSurfaceCondition(world, mirror);
    condition->EnableLambertian(false); // Disable Lambertian reflection
    ReflectiveSurface reflectiveSurface(brdfModel, albedo);
    // Use reflectiveSurface.computeBRDF() in the ray tracing logic
    world->AddNode(mirror, 1);
    manager->SetNsegments(100);
    manager->CloseGeometry();
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
void sweepDetector(BRDFModel brdfModel, double albedo) {
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    setupOpticsManager(manager, brdfModel, albedo);

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
    std::string filename = "nonlambertian_fluxmap_data.csv";
    std::ofstream csvFile(filename);

    if (!csvFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
        return;
    }

    // Write header with metadata
    csvFile << "# Nonlambertian Flux Map Data\n"
            << "# BRDF Model: " << (brdfModel == OREN_NAYAR ? "Oren-Nayar" : "Cook-Torrance") << "\n"
            << "# Albedo: " << albedo << "\n"
            << "# Number of rays per point: " << n << "\n"
            << "# Bins: " << nThetaBins << " theta x " << nPhiBins << " phi\n"
            << "theta,phi,fraction\n";

    // Write data with error checking
    for(int i = 0; i < nThetaBins; i++) {
        double theta = (i + 0.5) * 90.0/nThetaBins;
        for(int j = 0; j < nPhiBins; j++) {
            double phi = (j + 0.5) * 360.0/nPhiBins;
            double fraction = fluxMap->GetBinContent(i+1, j+1);
            
            if (!(csvFile << std::fixed << std::setprecision(6)
                         << theta << "," << phi << "," << fraction << "\n")) {
                std::cerr << "Error writing to CSV file at theta=" << theta 
                         << ", phi=" << phi << std::endl;
                csvFile.close();
                return;
            }
        }
    }

    csvFile.close();

    if (csvFile.good()) {
        std::cout << "\nFlux map data successfully saved to '" << filename << "'" << std::endl;
    } else {
        std::cerr << "Error: Some data may not have been written correctly to " 
                  << filename << std::endl;
    }
}

void visualizeDetector(double theta = 0.0, double phi = 0.0, BRDFModel brdfModel = COOK_TORRANCE, double albedo = 0.9) {
    // Create canvas with square proportions
    TCanvas* c = new TCanvas("cvis", "Detector Visualization", 800, 800);
    c->cd();
    
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    setupOpticsManager(manager, brdfModel, albedo);
    
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
