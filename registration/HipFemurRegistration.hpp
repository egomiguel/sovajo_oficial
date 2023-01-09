#ifndef HIP_FEMUR_REGISTRATION_H
#define HIP_FEMUR_REGISTRATION_H

#include "Registration.hpp"
#include "registration_export.h"

class REGISTRATION_EXPORT HipRegistrationFemur: public Registration
{
public:
    HipRegistrationFemur(const vtkSmartPointer<vtkPolyData> pImage, const PointTypeITK& pAnteriorFemoralNeckCT, const PointTypeITK& pAnteriorDistalTrochanterCT, const PointTypeITK& pLateralTrochanterCT, RegisterSide pSide);
    
    std::vector<RegistrationPointsHip> getRegistrationPointPosterolateral(double& pError) const;

    std::vector<RegistrationPointsHip> getRegistrationPointAnterolateral(double& pError) const;

    std::vector<RegistrationPointsHip> getRegistrationPointAnterior(double& pError) const;

    bool RegistrationLandmarks(const PointTypeITK& pAnteriorFemoralNeckCamera, const PointTypeITK& pAnteriorDistalTrochanterCamera, const PointTypeITK& pLateralTrochanterCamera, double& error);

    bool MakeRegistration(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pAnteriorFemoralNeckCamera, const PointTypeITK& pAnteriorDistalTrochanterCamera, const PointTypeITK& pLateralTrochanterCamera);

private:
    PointTypeITK mAnteriorFemoralNeck;
    PointTypeITK mAnteriorDistalTrochanter;
    PointTypeITK mLateralTrochanter;
    RegisterSide mSide;
    cv::Point3d mCenter;

    cv::Mat getTemplateAlignment(double& pError) const;
};

#endif