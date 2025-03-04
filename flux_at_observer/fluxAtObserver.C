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
        // Convert spherical to cartesian coordinates (relative to exit port at z=-100cm)
        double theta_rad = theta * M_PI / 180.0;
        double phi_rad = phi * M_PI / 180.0;
        
        // Calculate position relative to exit port
        x = radius * sin(theta_rad) * cos(phi_rad);
        y = radius * sin(theta_rad) * sin(phi_rad);
        z = -100*cm - radius * cos(theta_rad);
        
        // Normal vector always points to exit port center
        double dx = x;
        double dy = y;
        double dz = z + 100*cm;  // relative to exit port
        double mag = sqrt(dx*dx + dy*dy + dz*dz);
        
        nx = -dx/mag;  // Negative because we want to point towards the source
        ny = -dy/mag;
        nz = -dz/mag;
    }
    
    bool checkIntersection(const ARay* ray) {
        Double_t lastPoint[3];
        Double_t direction[3];
        ray->GetLastPoint(lastPoint);
        ray->GetDirection(direction); // Using GetDirection instead of GetLastStatus
        
        // Check if ray intersects detector plane
        double dx = lastPoint[0] - x;
        double dy = lastPoint[1] - y;
        double dz = lastPoint[2] - z;
        
        double dot = direction[0]*nx + direction[1]*ny + direction[2]*nz;
        if (dot >= 0) return false;  // Ray moving away from detector
        
        // Calculate intersection point
        double t = -(dx*nx + dy*ny + dz*nz) / dot;
        double ix = lastPoint[0] + direction[0] * t;
        double iy = lastPoint[1] + direction[1] * t;
        double iz = lastPoint[2] + direction[2] * t;
        
        // Check if intersection point is within detector bounds
        double local_x = (ix - x)*nx + (iy - y)*ny + (iz - z)*nz;
        double local_y = -(ix - x)*ny + (iy - y)*nx;
        
        return (fabs(local_x) <= width/2 && fabs(local_y) <= height/2);
    }

    TGeoVolume* CreateGeometry() const {
        // Create a thicker box to represent the detector
        TGeoBBox* detBox = new TGeoBBox("detBox", 
            width/2,    // half-width in x
            height/2,   // half-height in y
            2*cm        // Make it thicker for better visibility
        );
        TGeoVolume* detVol = new TGeoVolume("detector", detBox);
        detVol->SetLineColor(kBlue);
        detVol->SetLineWidth(2);
        detVol->SetFillColor(kCyan);          // Brighter color
        detVol->SetTransparency(20);          // More opaque
        return detVol;
    }
    
    void AddToGeometry(TGeoVolume* top) const {
        TGeoVolume* detVol = CreateGeometry();
        
        // Create rotation matrix to align detector normal
        // First, calculate rotation angles
        double phi_rot = atan2(ny, nx) * 180/M_PI;
        double theta_rot = acos(nz) * 180/M_PI;
        
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
    world->AddNode(mirror, 1);
    manager->SetNsegments(100);
    manager->CloseGeometry();
}

bool isRayPassingThroughExitPort(ARay* ray, double exitPortZ) {
    Double_t lastPoint[3];
    ray->GetLastPoint(lastPoint);
    return lastPoint[2] < exitPortZ;
}

void traceRays(AOpticsManager* manager, int n, double exitPortZ, Detector& detector) {
    for (int i = 0; i < n; ++i) {
        double x = -60*cm;
        ARay* ray = new ARay(0, 400*nm, x, 0*cm, -80*cm, 0, 5, 0, 0);
        manager->TraceNonSequential(*ray);
        TPolyLine3D* pol = ray->MakePolyLine3D();
        pol->SetLineWidth(1);
        pol->SetLineColor(2);
        pol->Draw();
        
        if (isRayPassingThroughExitPort(ray, exitPortZ) && 
            detector.checkIntersection(ray)) {
            detector.hitCount++;
        }

        delete ray;
    }
    std::cout << "Rays hitting detector: " << detector.hitCount << std::endl;
}

void sweepDetector() {
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    setupOpticsManager(manager);

    int n = 1000; // reduced number of rays per point for faster sweep
    double exitPortZ = -100*cm;
    
    // Create 2D histogram for the angular sweep
    const int nThetaBins = 45;  // 2-degree steps in theta
    const int nPhiBins = 30;    // 5-degree steps in phi
    TH2D* fluxMap = new TH2D("fluxMap", "Detector Flux Map;#theta (deg);#phi (deg)", 
                            nThetaBins, 0, 90, nPhiBins, 0, 360);
    
    // Create detector
    Detector detector;
    
    // Perform the sweep
    for(int i = 0; i < nThetaBins; i++) {
        double theta = (i + 0.5) * 90.0/nThetaBins;
        for(int j = 0; j < nPhiBins; j++) {
            double phi = (j + 0.5) * 360.0/nPhiBins;
            
            detector.hitCount = 0; // reset counter
            detector.setPosition(theta, phi, 100*cm);
            
            traceRays(manager, n, exitPortZ, detector);
            
            // Fill histogram with normalized counts
            fluxMap->SetBinContent(i+1, j+1, double(detector.hitCount)/n);
            
            // Progress update
            std::cout << "Completed theta = " << theta << "°, phi = " << phi << "°" << std::endl;
        }
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
    world->SetTransparency(80);  // Make world box very transparent
    
    // Add detector to geometry
    detector.AddToGeometry(world);
    
    // Draw the geometry with better visualization
    gGeoManager->CloseGeometry();
    world->Draw("ogl");
    
    // Trace rays
    int n = 50;
    double exitPortZ = -100*cm;
    traceRays(manager, n, exitPortZ, detector);
    
    //draw a line from the exit port to the detector
    TGeoTube* tube = new TGeoTube("tube", 0, 0.1*cm, 100*cm);
    TGeoVolume* line = new TGeoVolume("line", tube);
    line->SetLineColor(kRed);
    line->SetLineWidth(10);

    // Create and draw coordinate axes
    double axisLength = 50*cm;
    // X axis (red)
    TGeoTube* xTube = new TGeoTube("xTube", 0, 0.2*cm, axisLength/2);
    TGeoVolume* xAxis = new TGeoVolume("xAxis", xTube);
    xAxis->SetLineColor(kRed);
    gGeoManager->GetTopVolume()->AddNode(xAxis, 1, new TGeoRotation("xrot", 0, 90, 0));

    // Y axis (green)
    TGeoTube* yTube = new TGeoTube("yTube", 0, 0.2*cm, axisLength/2);
    TGeoVolume* yAxis = new TGeoVolume("yAxis", yTube);
    yAxis->SetLineColor(kGreen);
    gGeoManager->GetTopVolume()->AddNode(yAxis, 1, new TGeoRotation("yrot", -90, 0, 0));

    // Z axis (blue)
    TGeoTube* zTube = new TGeoTube("zTube", 0, 0.2*cm, axisLength/2);
    TGeoVolume* zAxis = new TGeoVolume("zAxis", zTube);
    zAxis->SetLineColor(kBlue);
    gGeoManager->GetTopVolume()->AddNode(zAxis, 1);

    c->Update();
    
    // Print detector information
    cout << "\nDetector Information:" << endl;
    cout << "Position (x,y,z): (" << detector.x/cm << ", " << detector.y/cm << ", " << detector.z/cm << ") cm" << endl;
    cout << "Normal vector (nx,ny,nz): (" << detector.nx << ", " << detector.ny << ", " << detector.nz << ")" << endl;
    cout << "Angular position: theta = " << theta << "°, phi = " << phi << "°" << endl;
}
