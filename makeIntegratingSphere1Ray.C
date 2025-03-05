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

const double cm = AOpticsManager::cm();
const double nm = AOpticsManager::nm();
const double thetaMax = 170.;     // governs size of exit port

void makeIntegratingSphere1Ray( )
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

    // Create and get canvas explicitly
    TCanvas* can = new TCanvas("can3D", "Integrating Sphere with Ray", 800, 800);
    can->cd();  // Make this the active canvas
    
    // Draw geometry first
    world->Draw("ogl");  // Use world instead of gGeoManager->GetTopVolume()
    
    // Create and trace ray
    ARay* ray = new ARay(0, 400*nm, -60*cm, 0*cm, -80*cm, 0, 5, 0, 0);
    manager->TraceNonSequential(*ray);
    TPolyLine3D* pol = ray->MakePolyLine3D();
    pol->SetLineWidth(1);
    pol->SetLineColor(2);
    
    // Draw ray in same view
    pol->Draw("same");
    
    // Update canvas
    can->Update();
    
    // Optional: delete ray;
}
