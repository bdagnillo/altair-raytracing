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

// Fix traceRays function to update detector's hit count
int traceRays(AOpticsManager* manager, int n, double exitPortZ, Detector& detector, bool drawRays = false) {
    std::vector<TPolyLine3D*> rays;  // Store rays for visualization
    int hitCount = 0;  // Count of rays hitting the detector
    
    // Draw a marker at the ray source position
    if (drawRays) {
        // Source position coordinates
        double sourceX = -60*cm;
        double sourceY = 0*cm;
        double sourceZ = -80*cm;
        
        // Create a 3D marker at the ray source
        TPolyMarker3D* sourceMarker = new TPolyMarker3D(1);
        sourceMarker->SetPoint(0, sourceX, sourceY, sourceZ);
        sourceMarker->SetMarkerColor(kMagenta);
        sourceMarker->SetMarkerStyle(20);  // Filled circle
        sourceMarker->SetMarkerSize(2.5);
        sourceMarker->Draw();
        
        // Print information about the source position
        std::cout << "Ray source position: (" << sourceX/cm << ", " 
                  << sourceY/cm << ", " << sourceZ/cm << ") cm" << std::endl;
    }
    
    for (int i = 0; i < n; ++i) {
        double x = -60*cm;
        double y = 0*cm;
        double z = -80*cm;
        // Direction
        double dirX = 5;
        double dirY = 2;
        double dirZ = 0;
        ARay* ray = new ARay(0, 660*nm, x, y, z, 0, dirX, dirY, dirZ);
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

    int n = 50000; // reduced number of rays per point for faster sweep
    double exitPortZ = -100*cm;
    
    // Create 2D histogram for the angular sweep
    const int nThetaBins = 180;
    const int nPhiBins = 90;
    TH2D* fluxMap = new TH2D("fluxMap", "Detector Flux Map;#theta (deg);#phi (deg)", 
                            nThetaBins, 0, 90, nPhiBins, 0, 360);
    
    // Create detector
    Detector detector;
    
    // ----- SETUP CSV FILE BEFORE THE SWEEP BEGINS -----
    // Save flux map data to CSV with improved format
    const char* saveFolder = "results";
    
    // Create directory for results
    std::string mkdirCmd = "mkdir -p \"" + std::string(saveFolder) + "\"";
    int result = system(mkdirCmd.c_str());
    if (result != 0) {
        std::cerr << "Warning: Could not create directory: " << saveFolder << std::endl;
    } else {
        std::cout << "Using directory: " << saveFolder << std::endl;
    }
    
    // Generate filename based on parameters
    std::string filename = "fluxmap_data_" + std::to_string(n) + "rays_" + 
                          std::to_string(nThetaBins * nPhiBins) + "points.csv";
    std::string fullPath = std::string(saveFolder) + "/" + filename;
    
    // Function to get a unique filename if the target file already exists
    auto getUniqueFilename = [](const std::string& basePath) -> std::string {
        // First check if file exists
        FILE* file = fopen(basePath.c_str(), "r");
        if (!file) {
            // File doesn't exist, so we can use this name
            return basePath;
        }
        fclose(file);
        
        // If file exists, find a unique name
        std::string directory;
        std::string filename;
        
        // Parse the path
        size_t lastSlash = basePath.find_last_of("/\\");
        if (lastSlash != std::string::npos) {
            directory = basePath.substr(0, lastSlash + 1);
            filename = basePath.substr(lastSlash + 1);
        } else {
            directory = "";
            filename = basePath;
        }
        
        // Split filename into stem and extension
        size_t lastDot = filename.find_last_of('.');
        std::string stem;
        std::string extension;
        if (lastDot != std::string::npos) {
            stem = filename.substr(0, lastDot);
            extension = filename.substr(lastDot);
        } else {
            stem = filename;
            extension = "";
        }
        
        // Find unique filename by incrementing counter
        int counter = 1;
        std::string newPath;
        do {
            newPath = directory + stem + "_" + std::to_string(counter) + extension;
            counter++;
            file = fopen(newPath.c_str(), "r");
            if (file) {
                fclose(file);
            } else {
                // Found a name that doesn't exist
                break;
            }
        } while (true);
        
        return newPath;
    };
    
    // Get unique filename if already exists
    fullPath = getUniqueFilename(fullPath);
    
    // Open CSV file with error checking
    std::ofstream csvFile(fullPath);
    if (!csvFile.is_open()) {
        std::cerr << "Error: Could not open file " << fullPath << " for writing." << std::endl;
        return;
    }
    
    // Write metadata as comments in CSV
    time_t now = time(nullptr);
    struct tm* timeinfo = localtime(&now);
    char timeBuffer[80];
    strftime(timeBuffer, sizeof(timeBuffer), "%Y-%m-%d %H:%M:%S", timeinfo);
    
    csvFile << "# Flux Map Data - Generated: " << timeBuffer << std::endl;
    csvFile << "# Number of rays per position: " << n << std::endl;
    csvFile << "# Detector dimensions: 10cm x 10cm" << std::endl;
    csvFile << "# Sphere inner radius: 100.1cm" << std::endl;
    csvFile << "# Sphere outer radius: 101cm" << std::endl;
    csvFile << "# Exit port angle: " << thetaMax << " degrees" << std::endl;
    csvFile << "# Theta bins: " << nThetaBins << std::endl;
    csvFile << "# Phi bins: " << nPhiBins << std::endl;
    csvFile << "# y direction: 2" << std::endl;
    csvFile << "theta,phi,fraction" << std::endl;  // Header
    
    std::cout << "\nStarting detector sweep..." << std::endl;
    std::cout << "Format: theta(°), phi(°): hits/total = fraction" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "CSV file created: " << fullPath << std::endl;
    
    // ----- PERFORM THE SWEEP AND WRITE DATA AS WE GO -----
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
            
            // Write data to CSV file immediately
            csvFile << std::fixed << std::setprecision(6)
                   << theta << "," << phi << "," << fraction << std::endl;
            // Flush to ensure data is written even if program crashes
            csvFile.flush();
        }
        // Add separator between theta values
        std::cout << "----------------------------------------" << std::endl;
    }
    
    // Add completion timestamp
    csvFile << "# Sweep completed at: " << timeBuffer << std::endl;
    csvFile.close();
    
    std::cout << "\nFlux map data saved to '" << fullPath << "'" << std::endl;
    
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
}

void visualizeDetector(double theta = 45.0, double phi = 0.0) {
    // Create canvas with square proportions
    TCanvas* c = new TCanvas("cvis", "Detector Visualization", 800, 800);
    c->cd();
    
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    
    // Setup basic geometry without closing it yet
    manager->SetLimit(10000);
    TGeoBBox* box = new TGeoBBox("box", 200*cm, 200*cm, 200*cm);
    AOpticalComponent* world = new AOpticalComponent("world", box);
    manager->SetTopVolume(world);
    TGeoSphere* sphere = new TGeoSphere("sphereWithExitPort", 100.1*cm, 101*cm, 0., thetaMax);
    AMirror* mirror = new AMirror("mirror", sphere);
    ABorderSurfaceCondition* condition = new ABorderSurfaceCondition(world, mirror);
    condition->EnableLambertian(true);
    world->AddNode(mirror, 1);
    manager->SetNsegments(100);
    
    // Create and position detector
    Detector detector(20*cm, 20*cm);
    detector.setPosition(theta, phi, 100*cm);
    
    // Add detector to geometry
    detector.AddToGeometry(world);
    
    // Now close the geometry after all components are added
    manager->CloseGeometry();
    
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
    int n = 100;
    double exitPortZ = -100*cm;
    int hitCount = traceRays(manager, n, exitPortZ, detector, true);
    
    // Final viewer updates - AVOID DrawGuides() which causes segfault
    if (glv) {
        // glv->DrawGuides(); - Remove this line that caused the crash
        glv->UpdateScene();
        glv->ResetCamerasAfterNextUpdate();
    }
    
    c->Update();
    
    // Print detector information
    cout << "\nDetector Information:" << endl;
    cout << "Position (x,y,z): (" << detector.x/cm << ", " << detector.y/cm << ", " << detector.z/cm << ") cm" << endl;
    cout << "Normal vector (nx,ny,nz): (" << detector.nx << ", " << detector.ny << ", " << detector.nz << ")" << endl;
    cout << "Angular position: theta = " << theta << "°, phi = " << phi << "°" << endl;
    cout << "Number of rays hitting the detector: " << hitCount << "/" << n << endl;
}

