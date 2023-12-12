#include "HipPelvisCupImplantMatch.hpp"
#include "ImplantsException.hpp"
#include "ImplantTools.hpp"
#include <itkVersorRigid3DTransform.h>

using namespace THA::IMPLANTS;

HipPelvisCupImplantMatch::HipPelvisCupImplantMatch()
{
    isInit = false;
}

void HipPelvisCupImplantMatch::init(const HipPelvis& pPelvis, const HipPelvisCupImplant& pImplant, const Point& pHipCenterOfRotation)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_PELVIS_CUT_IMPLANT_MATCH;
    }

    this->mPelvis = pPelvis;
    this->mImplant = pImplant;
    this->mHipCenterOfRotation = pHipCenterOfRotation;
    isInit = true;
}

//itk::Matrix< double, 3, 3 > HipPelvisCupImplantMatch::GetRotationMatrix() const
//{
//    itk::Matrix< double, 3, 3 > rotation;
//
//    rotation[0][0] = rotationMatrix.at <double>(0, 0);
//    rotation[0][1] = rotationMatrix.at <double>(0, 1);
//    rotation[0][2] = rotationMatrix.at <double>(0, 2);
//
//    rotation[1][0] = rotationMatrix.at <double>(1, 0);
//    rotation[1][1] = rotationMatrix.at <double>(1, 1);
//    rotation[1][2] = rotationMatrix.at <double>(1, 2);
//
//    rotation[2][0] = rotationMatrix.at <double>(2, 0);
//    rotation[2][1] = rotationMatrix.at <double>(2, 1);
//    rotation[2][2] = rotationMatrix.at <double>(2, 2);
//
//    return rotation;
//}
//
//itk::Vector< double, 3 > HipPelvisCupImplantMatch::GetTranslationMatrix() const
//{
//    itk::Vector< double, 3 > translation;
//
//    translation[0] = translationMatrix.at <double>(0, 0);
//    translation[1] = translationMatrix.at <double>(1, 0);
//    translation[2] = translationMatrix.at <double>(2, 0);
//
//    return translation;
//}

itk::Rigid3DTransform<>::Pointer HipPelvisCupImplantMatch::getTransform(double pAbductionAngle, double pAnteversionAngle, double pShifSuperior, double pShifLateral, double pShiftAnterior) const
{
    std::vector<cv::Point3d> implantVectors;
    std::vector<cv::Point3d> pelvisVectors;

    Point implantX = mImplant.getVectorX();
    Point implantZ = mImplant.getVectorZ();
    Point implantY = implantX.cross(implantZ);
    implantY.normalice();

    implantVectors.push_back(implantX.ToCVPoint());
    implantVectors.push_back(implantZ.ToCVPoint());
    implantVectors.push_back(implantY.ToCVPoint());

    auto abdAnt = mPelvis.getAbductionAnteversionVectorsZX(mHipCenterOfRotation, pAbductionAngle, pAnteversionAngle);

    Point pelvisX = abdAnt.second;
    Point pelvisZ = abdAnt.first;
    Point pelvisY = pelvisX.cross(pelvisZ);
    pelvisY.normalice();

    pelvisVectors.push_back(pelvisX.ToCVPoint());
    pelvisVectors.push_back(pelvisZ.ToCVPoint());
    pelvisVectors.push_back(pelvisY.ToCVPoint());

    cv::Mat implantMatrix = cv::Mat(implantVectors.size(), 3, CV_64F, implantVectors.data());
    cv::Mat pelvisMatrix = cv::Mat(pelvisVectors.size(), 3, CV_64F, pelvisVectors.data());

    cv::Mat inverse = (implantMatrix.t()).inv();
    cv::Mat rotationMatrix = (pelvisMatrix.t()) * inverse;
	cv::Mat translationMatrix = mHipCenterOfRotation.ToMatPoint() - (rotationMatrix * mImplant.getCenterOfRotationImplant().ToMatPoint());

	
	cv::Mat cupCenterMat = (rotationMatrix * mImplant.getCenterOfRotationImplant().ToMatPoint()) + translationMatrix;
	Point cupCenter = Point(cupCenterMat);
	
	Point newCenterCup = mHipCenterOfRotation + pShifSuperior * mPelvis.getPelvisVectorInfSup();
	newCenterCup = newCenterCup + pShifLateral * mPelvis.getPelvisVectorLateralASIS();
	newCenterCup = newCenterCup + pShiftAnterior * mPelvis.getPelvisVectorAP();

	Point diff = newCenterCup - cupCenter;
	translationMatrix = translationMatrix + diff.ToMatPoint();
	
	itk::Matrix< double, 3, 3 > rotation;
	itk::Vector< double, 3 > translation;

	translation[0] = translationMatrix.at <double>(0, 0);
	translation[1] = translationMatrix.at <double>(1, 0);
	translation[2] = translationMatrix.at <double>(2, 0);

	rotation[0][0] = rotationMatrix.at <double>(0, 0);
	rotation[0][1] = rotationMatrix.at <double>(0, 1);
	rotation[0][2] = rotationMatrix.at <double>(0, 2);

	rotation[1][0] = rotationMatrix.at <double>(1, 0);
	rotation[1][1] = rotationMatrix.at <double>(1, 1);
	rotation[1][2] = rotationMatrix.at <double>(1, 2);

	rotation[2][0] = rotationMatrix.at <double>(2, 0);
	rotation[2][1] = rotationMatrix.at <double>(2, 1);
	rotation[2][2] = rotationMatrix.at <double>(2, 2);

	itk::Rigid3DTransform<double>::Pointer transform = itk::VersorRigid3DTransform<double>::New();

	transform->SetMatrix(rotation);
	transform->SetOffset(translation);

	return transform;
}

void HipPelvisCupImplantMatch::GetRobotTransform(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut) const
{
	Point cupCenter = mImplant.getCenterOfRotationImplant();
	cupCenter = ImplantTools::TransformPoint(cupCenter, pTransformIn);

	Plane base = mPelvis.getCoronalPlaneAPP();
	base.movePlane(cupCenter);

	Plane sagital;
	sagital.init(mPelvis.getPelvisVectorASIS(), mPelvis.getPubicJoin());
	sagital.reverseByPoint(mHipCenterOfRotation);

    Point vectorX = base.getProjectionVector(sagital.getNormalVector());
    Point vectorY = base.getProjectionVector(mPelvis.getPelvisVectorInfSup());

    cv::Mat rotationX = ImplantTools::GetRotateX(vectorX);

    cv::Mat tempRotationY = rotationX * vectorY.ToMatPoint();

    cv::Mat rotationY = ImplantTools::GetRotateY(Point(tempRotationY));

    cv::Mat myRotation = rotationY * rotationX;

    itk::Vector< double, 3 > translate;

    Point tTemp = myRotation * (cupCenter.ToMatPoint());

    translate[0] = -tTemp.x;
    translate[1] = -tTemp.y;
    translate[2] = -tTemp.z;

    itk::Matrix< double, 3, 3 > rotationITK;
    rotationITK[0][0] = myRotation.at<double>(0, 0);
    rotationITK[0][1] = myRotation.at<double>(0, 1);
    rotationITK[0][2] = myRotation.at<double>(0, 2);

    rotationITK[1][0] = myRotation.at<double>(1, 0);
    rotationITK[1][1] = myRotation.at<double>(1, 1);
    rotationITK[1][2] = myRotation.at<double>(1, 2);

    rotationITK[2][0] = myRotation.at<double>(2, 0);
    rotationITK[2][1] = myRotation.at<double>(2, 1);
    rotationITK[2][2] = myRotation.at<double>(2, 2);

    pTransformOut->SetMatrix(rotationITK);
    pTransformOut->SetOffset(translate);
}