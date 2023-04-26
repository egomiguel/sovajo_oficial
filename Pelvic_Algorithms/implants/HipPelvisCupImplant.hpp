#ifndef HIP_PELVIS_CUP_IMPLANT_H
#define HIP_PELVIS_CUP_IMPLANT_H

#include "implants_export.h"
#include "Plane.hpp"


class IMPLANTS_EXPORT HipPelvisCupImplant
{
public:
    HipPelvisCupImplant();
    void init(const Point& pTopPoint, const Point& pBasePoint1, const Point& pBasePoint2, const Point& pBasePoint3);
    Point getVectorX() const;
    Point getVectorZ() const;
    Point getTopPoint() const;
    Point getCenterOfRotationImplant() const;
    Plane getBasePlane() const;

private:
    Point mTopPoint;
    Point mBasePoint1, mBasePoint2, mBasePoint3;
    bool isInit;
};




#endif