#include "TCanvas.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoSphere.h"
#include "TGeoTube.h"
#include "TH2.h"
#include "TThread.h"
#include "TH1D.h"
#include "TMath.h"  // Add TMath header
#include <cmath>

#include "AFocalSurface.h"
#include "AGeoUtil.h"
#include "AMirror.h"
#include "AObscuration.h"
#include "AOpticsManager.h"
#include "ABorderSurfaceCondition.h"
#include "ARayShooter.h"
#include <iostream>

const double cm = AOpticsManager::cm();
const double nm = AOpticsManager::nm();
const double thetaMax = 170.;     // governs size of exit port

//plot angular distribution at exit port

void distributionSphereDetectorSweep( )
{
//    TThread::Initialize();
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    manager->SetLimit(10000);
//    manager->SetLimit(1);
//    TGeoBBox* box = new TGeoBBox("box", 100*cm, 100*cm, 100*cm);
    TGeoBBox* box = new TGeoBBox("box", 200*cm, 200*cm, 200*cm);
    AOpticalComponent* world = new AOpticalComponent("world", box);
    manager->SetTopVolume(world);
//    TGeoSphere* sphere = new TGeoSphere("sphere", 99*cm, 100*cm);
    TGeoSphere* sphere = new TGeoSphere("sphereWithExitPort", 100.1*cm, 101*cm, 0., thetaMax);
    AMirror* mirror = new AMirror("mirror", sphere);
    ABorderSurfaceCondition* condition = new ABorderSurfaceCondition(world, mirror);
    condition->EnableLambertian(true);
    world->AddNode(mirror, 1);
    manager->SetNsegments(100);
    manager->CloseGeometry();
    gGeoManager->GetTopVolume()->Draw("ogl");

    // Create histogram with more bins for better resolution
    TH1D* hAngularDist = new TH1D("hAngularDist", 
                                 "Angular Distribution of Exiting Rays;Angle from normal (degrees);Count", 
                                 180, -90, 90);  // Increased bin count
    
    // Add debug histograms for raw direction components
    TH2D* hDirectionsXZ = new TH2D("hDirectionsXZ", "Ray Direction Components X-Z;X;Z", 100, -1, 1, 100, -1, 1);
    TH2D* hDirectionsYZ = new TH2D("hDirectionsYZ", "Ray Direction Components Y-Z;Y;Z", 100, -1, 1, 100, -1, 1);
    TH1D* hDirectionZ = new TH1D("hDirectionZ", "Z Direction Component;Z;Count", 100, -1, 1);

    int n = 10000; // specify the number of rays
    int fluxCount = 0; // counter for rays passing through the exit port
    double exitPortZ = -100*cm; // z position of the exit port

    for (int i = 0; i < n; ++i) {
        double x = -60*cm;
        ARay* ray = new ARay(0, 400*nm, x, 0*cm, -80*cm, 0, 5, 0, 0);
        manager->TraceNonSequential(*ray);
        TPolyLine3D* pol = ray->MakePolyLine3D();
        pol->SetLineWidth(1);
        pol->SetLineColor(2);
        pol->Draw();

        // Get ray's last point and direction
        Double_t lastPoint[3] = {0, 0, 0};
        Double_t direction[3] = {0, 0, 0};
        ray->GetLastPoint(lastPoint);
        ray->GetDirection(direction);

        if (lastPoint[2] < exitPortZ) {
            fluxCount++;
            
            // Normalize direction vector first
            double norm = std::sqrt(direction[0]*direction[0] + 
                                  direction[1]*direction[1] + 
                                  direction[2]*direction[2]);
            
            double dx = direction[0]/norm;
            double dy = direction[1]/norm;
            double dz = direction[2]/norm;
            
            // Fill debug histograms
            hDirectionsXZ->Fill(dx, dz);
            hDirectionsYZ->Fill(dy, dz);
            hDirectionZ->Fill(dz);

            // Calculate angle from z axis
            double theta = TMath::Sign(std::acos(dz) * 180.0/M_PI, dx);  // angle from z axis
            
            // Fill histogram with proper weighting
            if (std::isfinite(theta)) {
                hAngularDist->Fill(theta, 1.0);
            }
        }

        delete ray;
    }

    std::cout << "Flux of rays through the exit port: " << fluxCount << std::endl;

    // Draw distributions
    TCanvas* c2 = new TCanvas("c2", "Angular Distribution", 800, 800);
    c2->Divide(2, 2);
    
    c2->cd(1);
    hAngularDist->Draw();
    
    // Add cosine fit with proper normalization
    TF1* cosfit = new TF1("cosfit", "[0]*cos((x)*TMath::Pi()/180.)", -90, 90);
    double integral = hAngularDist->Integral("width");  // get area considering bin widths
    cosfit->SetParameter(0, integral);
    hAngularDist->Fit(cosfit, "R");
    
    c2->cd(2);
    hDirectionZ->Draw();
    
    c2->cd(3);
    hDirectionsXZ->Draw("colz");
    
    c2->cd(4);
    hDirectionsYZ->Draw("colz");
    
    c2->Update();
}
