#ifndef HIP_FEMUR_STEM_IMPLANT_H
#define HIP_FEMUR_STEM_IMPLANT_H

#include "implants_export.h"
#include "Plane.hpp"


class IMPLANTS_EXPORT HipFemurStemImplant
{
public:
    HipFemurStemImplant();
    void init(const Point& pTopPoint, const Point& pBasePoint, const Point& pHeadCenter);
    Point getVectorInfoSup() const;
    Point getVectorLatMed() const;
    Point getHeadCenter() const;

private:
    Point mTopPoint, mBasePoint, mHeadCenter;
    bool isInit;
};




#endif