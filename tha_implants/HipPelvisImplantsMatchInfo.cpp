#include "HipPelvisImplantsMatchInfo.hpp"
#include "ImplantsException.hpp"
#include "ImplantTools.hpp"
#include <itkVersorRigid3DTransform.h>

using namespace THA::IMPLANTS;

HipPelvisImplantsMatchInfo::HipPelvisImplantsMatchInfo(const HipPelvis& pPelvis, const HipPelvisCupImplant& pImplantCup, const Point& pHipCenterOfRotation, const itk::Rigid3DTransform<>::Pointer pImplantToBoneCupTransform)
{
	this->mPelvis = pPelvis;
	this->mImplantCup = pImplantCup;
	this->mHipCenterOfRotation = pHipCenterOfRotation;
	setCupTransform(pImplantToBoneCupTransform);
}

HipPelvisImplantsMatchInfo::~HipPelvisImplantsMatchInfo()
{
}

void HipPelvisImplantsMatchInfo::setCupTransform(const itk::Rigid3DTransform<>::Pointer pImplantToBoneStemTransform)
{
	mRotationCup = Rigid3DTransformToCVRotation(pImplantToBoneStemTransform);
	mTranslationCup = Rigid3DTransformToCVTranslation(pImplantToBoneStemTransform);
}

itk::Rigid3DTransform<>::Pointer HipPelvisImplantsMatchInfo::getITKCupTransform() const
{
	itk::Matrix< double, 3, 3 > rotation;
	itk::Vector< double, 3 > translate;

	rotation[0][0] = mRotationCup.at<double>(0, 0);
	rotation[0][1] = mRotationCup.at<double>(0, 1);
	rotation[0][2] = mRotationCup.at<double>(0, 2);

	rotation[1][0] = mRotationCup.at<double>(1, 0);
	rotation[1][1] = mRotationCup.at<double>(1, 1);
	rotation[1][2] = mRotationCup.at<double>(1, 2);

	rotation[2][0] = mRotationCup.at<double>(2, 0);
	rotation[2][1] = mRotationCup.at<double>(2, 1);
	rotation[2][2] = mRotationCup.at<double>(2, 2);

	translate[0] = mTranslationCup.at<double>(0, 0);
	translate[1] = mTranslationCup.at<double>(1, 0);
	translate[2] = mTranslationCup.at<double>(2, 0);

	itk::Rigid3DTransform<double>::Pointer transform = itk::VersorRigid3DTransform<double>::New();
	transform->SetMatrix(rotation);
	transform->SetOffset(translate);

	return transform;
}

cv::Mat HipPelvisImplantsMatchInfo::Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer transform) const
{
	itk::Matrix< double, 3, 3 > rotation = transform->GetMatrix();
	double* matrix = new double[9];

	matrix[0] = rotation[0][0];
	matrix[1] = rotation[0][1];
	matrix[2] = rotation[0][2];

	matrix[3] = rotation[1][0];
	matrix[4] = rotation[1][1];
	matrix[5] = rotation[1][2];

	matrix[6] = rotation[2][0];
	matrix[7] = rotation[2][1];
	matrix[8] = rotation[2][2];

	cv::Mat matrixCV(3, 3, CV_64FC1, matrix);
	return matrixCV;
}

cv::Mat HipPelvisImplantsMatchInfo::Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer transform) const
{
	itk::Vector< double, 3 > translate = transform->GetOffset();
	cv::Mat result(3, 1, CV_64FC1);

	result.at <double>(0, 0) = translate[0];
	result.at <double>(1, 0) = translate[1];
	result.at <double>(2, 0) = translate[2];
	return result;
}

double HipPelvisImplantsMatchInfo::getCupInclination() const
{
	Point rotationAxis = mPelvis.getPelvisVectorAP();
	Point referenceVector = mPelvis.getPelvisVectorInfSup();

	cv::Mat implantVectorMat = mRotationCup * mImplantCup.getVectorZ().ToMatPoint();
	Point implantVector = Point(implantVectorMat);

	Plane coronal;
	coronal.init(rotationAxis, mHipCenterOfRotation);

	Point vectorImplantProj = coronal.getProjectionVector(implantVector);
	Point vectorBoneProj = coronal.getProjectionVector(referenceVector);
	vectorImplantProj.normalice();
	vectorBoneProj.normalice();

	double angle = ImplantTools::getAngleBetweenVectorsDegree(vectorImplantProj, vectorBoneProj);

	Plane sagital;
	Point ref = mPelvis.getPubicJoin();
	sagital.init(mPelvis.getPelvisVectorLateralASIS(), ref);
	sagital.reverseByPoint(mHipCenterOfRotation, false);

	ref = ref + 1000. * vectorImplantProj;

	if (sagital.eval(ref) < 0)
	{
		return -angle;
	}
	else
	{
		return angle;
	}
}

double HipPelvisImplantsMatchInfo::getCupAntversion() const
{
	cv::Mat implantVectorMat = mRotationCup * mImplantCup.getVectorZ().ToMatPoint();
	Point implantVector = Point(implantVectorMat);
	implantVector.normalice();

	Point refPoint = mPelvis.getCoronalPlaneAPP().getProjectionPoint(mPelvis.getPubicJoin());

	Point pelvisVectorAP = mPelvis.getPelvisVectorAP();

	/*
		The plane formed by the normal of the implant and the AP vector of the pelvis is used.
		The anteversion rotation occurs in that plane, therefore in that plane the angle is sought
	*/

	Point interceptionVector = mPelvis.getCoronalPlaneAPP().getInterceptionPlaneVector(implantVector, pelvisVectorAP);
	interceptionVector.normalice();

	double angle1 = ImplantTools::getAngleBetweenVectorsDegree(interceptionVector, implantVector);
	double angle2 = ImplantTools::getAngleBetweenVectorsDegree(-interceptionVector, implantVector);
	double angle = std::min(angle1, angle2);

	refPoint = refPoint + 1000. * implantVector;

	if (mPelvis.getCoronalPlaneAPP().eval(refPoint) < 0)
	{
		return angle;
	}
	else
	{
		return -angle;
	}
}

double HipPelvisImplantsMatchInfo::getCupShiftSuperior() const
{
	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;

	Line superiorLine(mPelvis.getPelvisVectorInfSup(), mHipCenterOfRotation);

	Plane axial;
	axial.init(mPelvis.getPelvisVectorInfSup(), mHipCenterOfRotation);

	Point cupCenter = Point(cupCenterMat);
	cupCenter = superiorLine.getProjectPoint(cupCenter);
	double distance = ImplantTools::getDistanceBetweenPoints(cupCenter, mHipCenterOfRotation);

	if (axial.eval(cupCenter) >= 0)
	{
		return distance;
	}
	else
	{
		return -distance;
	}
}

double HipPelvisImplantsMatchInfo::getCupShiftLateral() const
{
	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;

	Line lateralLine(mPelvis.getPelvisVectorLateralASIS(), mHipCenterOfRotation);

	Plane sagital;
	sagital.init(mPelvis.getPelvisVectorLateralASIS(), mHipCenterOfRotation);

	Point cupCenter = Point(cupCenterMat);
	cupCenter = lateralLine.getProjectPoint(cupCenter);
	double distance = ImplantTools::getDistanceBetweenPoints(cupCenter, mHipCenterOfRotation);

	if (sagital.eval(cupCenter) >= 0)
	{
		return distance;
	}
	else
	{
		return -distance;
	}
}

double HipPelvisImplantsMatchInfo::getCupShiftAnterior() const
{
	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;

	Line anteriorLine(mPelvis.getPelvisVectorAP(), mHipCenterOfRotation);

	Plane coronal;
	coronal.init(mPelvis.getPelvisVectorAP(), mHipCenterOfRotation);

	Point cupCenter = Point(cupCenterMat);
	cupCenter = anteriorLine.getProjectPoint(cupCenter);
	double distance = ImplantTools::getDistanceBetweenPoints(cupCenter, mHipCenterOfRotation);

	if (coronal.eval(cupCenter) >= 0)
	{
		return distance;
	}
	else
	{
		return -distance;
	}
}


Plane HipPelvisImplantsMatchInfo::TransformPlane(const Plane& plane, const cv::Mat& rotation, const cv::Mat& translation) const
{
	cv::Mat transformNormalVector = rotation * plane.getNormalVectorMat();
	cv::Mat transformPoint = (rotation * plane.getPointMat()) + translation;
	Plane transformPlane;
	transformPlane.init(Point(transformNormalVector), Point(transformPoint));
	return transformPlane;
}