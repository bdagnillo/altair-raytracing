#ifndef REFLECTIVE_SURFACE_H
#define REFLECTIVE_SURFACE_H

#include "TVector3.h"
#include "ARay.h"

enum BRDFModel {
    OREN_NAYAR,
    COOK_TORRANCE
};

class ReflectiveSurface {
public:
    ReflectiveSurface(BRDFModel model, Double_t albedo);
    
    TVector3 sampleHemisphereCosine(Double_t u1, Double_t u2, const TVector3& normal) const;
    Double_t computeBRDF(const TVector3& incoming, const TVector3& outgoing, const TVector3& normal) const;
    
    static TVector3 getRayDirection(const ARay* ray);
    static void setRayDirection(ARay& ray, const TVector3& direction);

private:
    BRDFModel fModel;
    Double_t fAlbedo;
    
    Double_t orenNayarBRDF(const TVector3& incoming, const TVector3& outgoing, const TVector3& normal) const;
    Double_t cookTorranceBRDF(const TVector3& incoming, const TVector3& outgoing, const TVector3& normal) const;
};

#endif
