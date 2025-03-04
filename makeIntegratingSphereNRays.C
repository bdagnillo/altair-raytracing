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
    gGeoManager->GetTopVolume()->Draw("ogl");

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

        // check if the ray passes through the exit port
        Double_t lastPoint[3] = {0, 0, 0};
        ray->GetLastPoint(lastPoint);
        // std::cout << "Ray " << i << " last point: (" << lastPoint[0] << ", " << lastPoint[1] << ", " << lastPoint[2] << ")" << std::endl;
        if (lastPoint[2] < exitPortZ) {
            fluxCount++;
        }

        delete ray;
    }

    std::cout << "Flux of rays through the exit port: " << fluxCount << std::endl;
}
