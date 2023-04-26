#ifndef HIP_FEMUR_STEM_IMPLANT_MATCH_H
#define HIP_FEMUR_STEM_IMPLANT_MATCH_H

#include "HipFemurStemImplant.hpp"
#include "HipPelvis.hpp"
#include <itkRigid3DTransform.h>
#include "implants_export.h"


class IMPLANTS_EXPORT HipFemurStemImplantMatch
{
public:
    HipFemurStemImplantMatch();

    void init(const HipPelvis& pPelvis, const HipFemurStemImplant& pImplant, const Point& pHipCenterOfRotation);

    itk::Matrix< double, 3, 3 > GetRotationMatrix() const;

    itk::Vector< double, 3 > GetTranslationMatrix() const;

private:
    HipPelvis mPelvis;
    HipFemurStemImplant mImplant;
    Point mHipCenterOfRotation;
    cv::Mat rotationMatrix;
    cv::Mat translationMatrix;

    void getRigidTransform();

    bool isInit;
};




#endif