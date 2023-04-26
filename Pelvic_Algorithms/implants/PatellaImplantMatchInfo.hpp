#ifndef PATELLA_IMPLANT_MATCH_INFO_H
#define PATELLA_IMPLANT_MATCH_INFO_H

#include "PatellaImplant.hpp"
#include "Knee.hpp"
#include <itkRigid3DTransform.h>
#include "implants_export.h"

class IMPLANTS_EXPORT PatellaImplantMatchInfo
{
public:
    PatellaImplantMatchInfo();

	void init(const PatellaImplant& implant, const Knee& knee, const itk::Rigid3DTransform<>::Pointer pImplantToBonePatellaTransform);

    void setPatellaTransform(const itk::Rigid3DTransform<>::Pointer pImplantToBonePatellaTransform);

    double getCutThickness() const;

    double getAngleML() const;

    double getAngleSI() const;

private:
    PatellaImplant implant;
	Knee knee;
	bool isInit;

    cv::Mat patellaRotation, patellaTranslation;

    cv::Mat Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer transform) const;

    cv::Mat Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer transform) const;

    Plane TransformPlane(const Plane& plane, const cv::Mat& rotation, const cv::Mat& translation) const;

};


#endif
