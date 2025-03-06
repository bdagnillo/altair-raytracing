#include "CustomBRDFMirror.h"
#include "TRandom.h"

ClassImp(CustomBRDFMirror)

CustomBRDFMirror::CustomBRDFMirror() : AMirror(), fSurface(COOK_TORRANCE, 0.9) {}

CustomBRDFMirror::CustomBRDFMirror(const char* name, TGeoShape* shape, BRDFModel model, Double_t albedo)
    : AMirror(name, shape), fSurface(model, albedo) {}

CustomBRDFMirror::~CustomBRDFMirror() {}

Bool_t CustomBRDFMirror::Reflect(ARay& ray, TVector3& normal) {
    TVector3 incoming = ReflectiveSurface::getRayDirection(&ray);
    
    // Sample new direction using BRDF
    Double_t u1 = gRandom->Uniform(0, 1);
    Double_t u2 = gRandom->Uniform(0, 1);
    TVector3 outgoing = fSurface.sampleHemisphereCosine(u1, u2, normal);
    
    // Set new direction
    ReflectiveSurface::setRayDirection(ray, outgoing);
    
    // Add point to continue ray path
    Double_t pos[3];
    ray.GetLastPoint(pos);
    ray.AddPoint(pos[0], pos[1], pos[2], 0.0);
    
    return kTRUE;
}
