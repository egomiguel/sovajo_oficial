#include "HipFemurStemImplantMatch.hpp"
#include "ImplantsException.hpp"


HipFemurStemImplantMatch::HipFemurStemImplantMatch()
{
    isInit = false;
}

void HipFemurStemImplantMatch::init(const HipPelvis& pPelvis, const HipFemurStemImplant& pImplant, const Point& pHipCenterOfRotation)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_FEMUR_STEM_IMPLANT_MATCH;
    }

    this->mPelvis = pPelvis;
    this->mImplant = pImplant;
    this->mHipCenterOfRotation = pHipCenterOfRotation;
    getRigidTransform();
    isInit = true;
}

itk::Matrix< double, 3, 3 > HipFemurStemImplantMatch::GetRotationMatrix() const
{
    itk::Matrix< double, 3, 3 > rotation;

    rotation[0][0] = rotationMatrix.at <double>(0, 0);
    rotation[0][1] = rotationMatrix.at <double>(0, 1);
    rotation[0][2] = rotationMatrix.at <double>(0, 2);

    rotation[1][0] = rotationMatrix.at <double>(1, 0);
    rotation[1][1] = rotationMatrix.at <double>(1, 1);
    rotation[1][2] = rotationMatrix.at <double>(1, 2);

    rotation[2][0] = rotationMatrix.at <double>(2, 0);
    rotation[2][1] = rotationMatrix.at <double>(2, 1);
    rotation[2][2] = rotationMatrix.at <double>(2, 2);

    return rotation;
}

itk::Vector< double, 3 > HipFemurStemImplantMatch::GetTranslationMatrix() const
{
    itk::Vector< double, 3 > translation;

    translation[0] = translationMatrix.at <double>(0, 0);
    translation[1] = translationMatrix.at <double>(1, 0);
    translation[2] = translationMatrix.at <double>(2, 0);

    return translation;
}

void HipFemurStemImplantMatch::getRigidTransform()
{
    std::vector<cv::Point3d> implantVectors;
    std::vector<cv::Point3d> femurVectors;

    Point implantX = mImplant.getVectorLatMed();
    Point implantZ = mImplant.getVectorInfoSup();
    Point implantY = implantX.cross(implantZ);
    implantY.normalice();

    implantVectors.push_back(implantX.ToCVPoint());
    implantVectors.push_back(implantZ.ToCVPoint());
    implantVectors.push_back(implantY.ToCVPoint());

    Point femurX = mPelvis.getFemurVectorLatMed(mHipCenterOfRotation);
    Point femurZ = mPelvis.getFemurVectorInfSup();
    Point femurY = femurX.cross(femurZ);
    femurY.normalice();

    femurVectors.push_back(femurX.ToCVPoint());
    femurVectors.push_back(femurZ.ToCVPoint());
    femurVectors.push_back(femurY.ToCVPoint());

    cv::Mat implantMatrix = cv::Mat(implantVectors.size(), 3, CV_64F, implantVectors.data());
    cv::Mat pelvisMatrix = cv::Mat(femurVectors.size(), 3, CV_64F, femurVectors.data());

    cv::Mat inverse = (implantMatrix.t()).inv();
    rotationMatrix = (pelvisMatrix.t()) * inverse;
    translationMatrix = mHipCenterOfRotation.ToMatPoint() - (rotationMatrix * mImplant.getHeadCenter().ToMatPoint());
}