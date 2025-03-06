#ifndef CUSTOM_MIRROR_H
#define CUSTOM_MIRROR_H

#include "AMirror.h"
#include "TVector3.h"
#include <functional>

///////////////////////////////////////////////////////////////////////////////
//
// CustomMirror
//
// A mirror class that allows custom reflection behavior through a callback
//
///////////////////////////////////////////////////////////////////////////////

class CustomMirror : public AMirror {
private:
    std::function<TVector3(const TVector3&, const TVector3&)> reflectionHandler;
    
public:
    CustomMirror(const char* name, TGeoShape* shape) : AMirror(name, shape), reflectionHandler(nullptr) {}
    virtual ~CustomMirror() {}
    
    void SetReflectionHandler(std::function<TVector3(const TVector3&, const TVector3&)> handler) {
        reflectionHandler = handler;
    }
    
    virtual bool Reflect(ARay& ray, const TVector3& normal) override {
        if (!reflectionHandler) {
            // Fall back to normal mirror behavior if no handler is set
            return AMirror::Reflect(ray, normal);
        }
        
        Double_t dir[3];
        ray.GetDirection(dir);
        TVector3 incoming(dir[0], dir[1], dir[2]);
        
        // Get new direction from handler
        TVector3 newDir = reflectionHandler(incoming, normal);
        
        // Update ray direction and add new point
        Double_t pos[3];
        ray.GetLastPoint(pos);
        ray.AddPoint(pos[0], pos[1], pos[2]);
        ray.SetDirection(newDir.X(), newDir.Y(), newDir.Z());
        
        return true;
    }
    
    ClassDef(CustomMirror, 1)
};

#endif // CUSTOM_MIRROR_H