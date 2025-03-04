#include "TCanvas.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoSphere.h"
#include "TGeoTube.h"
#include "TH2.h"
#include "TH1.h"
#include "TThread.h"
#include "TGraph.h"
#include "TMath.h"

#include "AFocalSurface.h"
#include "AGeoUtil.h"
#include "AMirror.h"
#include "AObscuration.h"
#include "AOpticsManager.h"
#include "ABorderSurfaceCondition.h"
#include "ARayShooter.h"
#include <iostream>
#include <cmath>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>

// flux through disk

// Forward declarations
bool isRayHittingDetector(ARay* ray);
void addDetectorDisk(AOpticalComponent* world, double theta, double phi, double diskRadius);

void sweepDetector(AOpticsManager* manager, AOpticalComponent* world, double diskRadius, 
                  int nRays, double dtheta, double thetaMax) {
    double dphi = 180;  // phi increment in degrees
    
    // Create output file for results
    std::ofstream outFile("detector_sweep3.txt");
    outFile << "Theta(deg)\tPhi(deg)\tHitFraction\n";
    
    // Create 2D histogram for theta-phi map
    TH2D* hSweepMap = new TH2D("hSweepMap", "Hit Fraction Map;Theta (deg);Phi (deg)", 
                               int(2*thetaMax/dtheta), -thetaMax, thetaMax,
                               int(360/dphi), 0, 360);

    // Sweep through theta and phi values
    for (double theta = -thetaMax; theta <= thetaMax; theta += dtheta) {
        for (double phi = 0; phi < 360; phi += dphi) {
            // Initial geometry setup
            if (theta == -thetaMax && phi == 0) {
                manager->SetNsegments(100);
                manager->CloseGeometry();
                gGeoManager->GetTopVolume()->Draw("ogl");
            } else {
                // Clear previous detector if exists
                TGeoNode* oldDetector = gGeoManager->GetTopVolume()->FindNode("detector_1");
                if (oldDetector) {
                    gGeoManager->GetTopVolume()->RemoveNode(oldDetector);
                }
            }
            
            // Add detector at new position
            addDetectorDisk(world, theta, phi, diskRadius);
            manager->CloseGeometry();

            int detectorHits = 0;
            
            // Perform ray tracing
            for (int i = 0; i < nRays; ++i) {
                double x = -60*AOpticsManager::cm();
                ARay* ray = new ARay(0, 400*AOpticsManager::nm(), x, 0, -80*AOpticsManager::cm(), 0, 5, 0, 0);
                manager->TraceNonSequential(*ray);
                
                if (isRayHittingDetector(ray)) {
                    detectorHits++;
                }
                
                delete ray;
            }

            double hitFraction = static_cast<double>(detectorHits) / nRays;
            std::cout << "Theta: " << theta << "° Phi: " << phi << "° Hit fraction: " << hitFraction << std::endl;
            outFile << theta << "\t" << phi << "\t" << hitFraction << "\n";
            hSweepMap->Fill(theta, phi, hitFraction);
        }
    }

    outFile.close();

    // Create plots of the results
    TCanvas* c3 = new TCanvas("c3", "Detector Sweep Results", 1200, 500);
    c3->Divide(2, 1);
    
    // Draw theta profile (averaged over phi)
    c3->cd(1);
    TH1D* hThetaProfile = hSweepMap->ProjectionX("hThetaProfile");
    hThetaProfile->Scale(1.0/int(360/dphi));  // Average over phi bins
    hThetaProfile->SetTitle("Average Hit Fraction vs Theta;Theta (degrees);Average Hit Fraction");
    hThetaProfile->Draw();
    
    // Draw 2D map
    c3->cd(2);
    hSweepMap->SetTitle("Hit Fraction Map");
    hSweepMap->Draw("colz");
    
    c3->Update();
}

void integratingSphereDetectorSweep()
{
    // Move constants inside function
    const double cm = AOpticsManager::cm();
    const double nm = AOpticsManager::nm();
    const double tmax = 170.;     // governs size of exit port

    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    manager->SetLimit(10000);
    TGeoBBox* box = new TGeoBBox("box", 200*cm, 200*cm, 200*cm);
    AOpticalComponent* world = new AOpticalComponent("world", box);
    manager->SetTopVolume(world);
    TGeoSphere* sphere = new TGeoSphere("sphereWithExitPort", 100.1*cm, 105*cm, 0., tmax);
    AMirror* mirror = new AMirror("mirror", sphere);
    ABorderSurfaceCondition* condition = new ABorderSurfaceCondition(world, mirror);
    condition->EnableLambertian(true);
    world->AddNode(mirror, 1);

    const int nRays = 100000;  // reduced number of rays per position due to more positions
    const double dtheta = 0.5;  // increased step size
    const double thetaMax = 45;
    double diskRadius = 5*cm;
    
    sweepDetector(manager, world, diskRadius, nRays, dtheta, thetaMax);
}

// Helper function implementations
bool isRayHittingDetector(ARay* ray) {
    const TObjArray* nodeHistory = ray->GetNodeHistory();
    for (int i = 0; i < nodeHistory->GetEntries(); ++i) {
        TGeoNode* node = (TGeoNode*)nodeHistory->At(i);
        if (TString(node->GetName()).Contains("detector")) {
            return true;
        }
    }
    return false;
}

void addDetectorDisk(AOpticalComponent* world, double theta, double phi, double diskRadius) {
    const double cm = AOpticsManager::cm();
    TGeoTube* detectorDisk = new TGeoTube("detectorDisk", 0, diskRadius, 0.1*cm);
    AFocalSurface* detector = new AFocalSurface("detector", detectorDisk);
    
    // Calculate position
    double r = 200*cm;
    double x = r * std::sin(theta * M_PI/180.0) * std::cos(phi * M_PI/180.0);
    double y = r * std::sin(theta * M_PI/180.0) * std::sin(phi * M_PI/180.0);
    double z = -r * std::cos(theta * M_PI/180.0); 
    
    // Calculate rotation angles to point towards exit port at (0,0,-100)
    double dx = 0 - x;
    double dy = 0 - y;
    double dz = -100*cm - z;
    
    // Calculate angles for rotation
    double rotTheta = -std::atan2(std::sqrt(dx*dx + dy*dy), dz) * 180.0/M_PI;
    double rotPhi = std::atan2(dy, dx) * 180.0/M_PI;
    
    // Create rotation matrix
    TGeoRotation* rot = new TGeoRotation("rot");
    rot->RotateZ(rotPhi);
    rot->RotateY(rotTheta);
    
    TGeoCombiTrans* transform = new TGeoCombiTrans("transform", x, y, z, rot);
    world->AddNode(detector, 1, transform);
}
