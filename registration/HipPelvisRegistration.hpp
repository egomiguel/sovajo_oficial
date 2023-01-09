#ifndef PELVIS_REGISTRATION_H
#define PELVIS_REGISTRATION_H

#include "Registration.hpp"
#include "registration_export.h"

class REGISTRATION_EXPORT PelvisRegistration: public Registration
{
public:
    PelvisRegistration(const vtkSmartPointer<vtkPolyData> pImage, const PointTypeITK& pAnteriorAcetabulum, const PointTypeITK& pPosteriorAcetabulum, const PointTypeITK& pSuperiorAcetabulum, RegisterSide pSide);
    
    std::vector<RegistrationPointsHip> getRegistrationPointPelvis(double& pError) const;

    bool RegistrationLandmarks(const PointTypeITK& pAnteriorAcetabulumCamera, const PointTypeITK& pPosteriorAcetabulumCamera, const PointTypeITK& pSuperiorAcetabulumCamera, double& error);

    bool MakeRegistration(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pAnteriorAcetabulumCamera, const PointTypeITK& pPosteriorAcetabulumCamera, const PointTypeITK& pSuperiorAcetabulumCamera);

private:
    PointTypeITK anteriorAcetabulum;
    PointTypeITK posteriorAcetabulum;
    PointTypeITK superiorAcetabulum;
    RegisterSide mSide;
    cv::Point3d mCenter;
    cv::Point3d mVectorNormal;
    double mRadius;

    /*std::vector<cv::Point3d> getAcetabularPoints();
    double getAxisReference(cv::Point3d& pVectorOutSideNormal, cv::Point3d& pVectorSuperiorRing, cv::Point3d& pCenter);*/
};

#endif