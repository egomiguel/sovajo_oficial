#ifndef HIP_PELVIS_CUP_IMPLANT_MATCH_H
#define HIP_PELVIS_CUP_IMPLANT_MATCH_H

#include "HipPelvisCupImplant.hpp"
#include "HipPelvis.hpp"
#include <itkRigid3DTransform.h>
#include "implants_export.h"


class IMPLANTS_EXPORT HipPelvisCupImplantMatch
{
public:
    HipPelvisCupImplantMatch();

    void init(const HipPelvis& pPelvis, const HipPelvisCupImplant& pImplant, const Point& pHipCenterOfRotation);

	itk::Rigid3DTransform<>::Pointer getTransform(double pAbductionAngle = 40, double pAnteversionAngle = 20);

	double getCupInclination(const itk::Rigid3DTransform<>::Pointer pTransform);

	double getCutVersion(const itk::Rigid3DTransform<>::Pointer pTransform);

    void GetRobotTransform(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut);

private:
    HipPelvis mPelvis;
    HipPelvisCupImplant mImplant;
    Point mHipCenterOfRotation;
    bool isInit;
};




#endif