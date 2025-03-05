#include "TCanvas.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoSphere.h"
#include "TGeoTube.h"
#include "TH2.h"
#include "TThread.h"

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

void makeIntegratingSphereNRays( )
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

    // Create a canvas and 3D viewer
    TCanvas* c1 = new TCanvas("c1", "Integrating Sphere with Rays", 800, 800);
    c1->cd();
    
    // Set up the 3D viewer more carefully
    TVirtualViewer3D* viewer = gPad->GetViewer3D("ogl");
    if (!viewer) {
        std::cerr << "Failed to create 3D viewer" << std::endl;
        return;
    }
    
    // Draw the geometry and wait for initialization
    gGeoManager->GetTopVolume()->Draw("ogl");
    gPad->Update();
    
    // Create rays after the viewer is ready
    int n = 1000; // reduce number of rays for testing
    int fluxCount = 0;
    double exitPortZ = -100*cm;
    
    // Create a TPolyMarker3D to store all ray paths
    TPolyLine3D** rays = new TPolyLine3D*[n];
    
    for (int i = 0; i < n; ++i) {
        double x = -60*cm;
        ARay* ray = new ARay(0, 400*nm, x, 0*cm, -80*cm, 0, 5, 0, 0);
        manager->TraceNonSequential(*ray);
        
        rays[i] = ray->MakePolyLine3D();
        rays[i]->SetLineWidth(1);
        rays[i]->SetLineColor(2);
        rays[i]->Draw();
        
        Double_t lastPoint[3] = {0, 0, 0};
        ray->GetLastPoint(lastPoint);
        if (lastPoint[2] < exitPortZ) {
            fluxCount++;
        }
        
        delete ray;
        
        // Update less frequently
        if (i % 500 == 0) {
            gPad->Modified();
            gPad->Update();
        }
    }
    
    // Final update
    gPad->Modified();
    gPad->Update();
    
    std::cout << "Flux of rays through the exit port: " << fluxCount << std::endl;
    
    // Clean up
    for (int i = 0; i < n; ++i) {
        delete rays[i];
    }
    delete[] rays;
}
