#ifndef CUSTOM_BRDF_MIRROR_H
#define CUSTOM_BRDF_MIRROR_H

#include "AMirror.h"
#include "TVector3.h"
#include "ReflectiveSurface.h"

class CustomBRDFMirror : public AMirror {
private:
    ReflectiveSurface fSurface;

public:
    CustomBRDFMirror();
    CustomBRDFMirror(const char* name, TGeoShape* shape, BRDFModel model, Double_t albedo);
    virtual ~CustomBRDFMirror();

    virtual Bool_t Reflect(ARay& ray, TVector3& normal);

    ClassDef(CustomBRDFMirror, 1)
};

#endif
