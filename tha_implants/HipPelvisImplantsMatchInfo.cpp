#include "HipPelvisImplantsMatchInfo.hpp"
#include "ImplantsException.hpp"
#include "ImplantTools.hpp"
#include <itkVersorRigid3DTransform.h>
#include "vtkSphereSource.h"
#include "vtkPlane.h"
#include "vtkPlaneCollection.h"
#include "vtkClipClosedSurface.h"

using namespace THA::IMPLANTS;

HipPelvisImplantsMatchInfo::HipPelvisImplantsMatchInfo(const HipPelvis& pPelvis, const Point& pHipCenterOfRotation, const HipPelvisCupImplant& pImplantCup,
	const HipFemurStemImplant& pImplantStem, const HipFemurStemHeadImplant& pImplantStemHead,
	const itk::Rigid3DTransform<>::Pointer pImplantToBoneCupTransform,
	const itk::Rigid3DTransform<>::Pointer pImplantToBoneStemTransform,
	const itk::Rigid3DTransform<>::Pointer pImplantHeadToStemTransform)
{
	this->mPelvis = pPelvis;
	this->mImplantCup = pImplantCup;
	this->mHipCenterOfRotation = pHipCenterOfRotation;
	this->mImplantStem = pImplantStem;
	this->mImplantStemHead = pImplantStemHead;
	setStemTransform(pImplantToBoneStemTransform);
	setCupTransform(pImplantToBoneCupTransform);
	setStemHeadTransform(pImplantHeadToStemTransform);
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
	implantVector.normalice();

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

double HipPelvisImplantsMatchInfo::getCupInclination(const itk::Rigid3DTransform<>::Pointer pImplantToBoneTransform) const
{
	Point rotationAxis = mPelvis.getPelvisVectorAP();
	Point referenceVector = mPelvis.getPelvisVectorInfSup();

	/*
		This vectorX represent the vector perpendicular to the base of the cup, 
		which should point towards the center of the hip and not outwards, 
		hence the change in sign. This transform is that of the robot when it uses the reamer.
	*/

	auto vectorX = pImplantToBoneTransform->GetMatrix().GetVnlMatrix().get_column(0);
	Point implantVectorZ = -Point(vectorX[0], vectorX[1], vectorX[2]);
	implantVectorZ.normalice();

	Plane coronal;
	coronal.init(rotationAxis, mHipCenterOfRotation);

	Point vectorImplantProj = coronal.getProjectionVector(implantVectorZ);
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

double HipPelvisImplantsMatchInfo::getCupAntversion(const itk::Rigid3DTransform<>::Pointer pImplantToBoneTransform) const
{
	/*
		This vectorX represent the vector perpendicular to the base of the cup,
		which should point towards the center of the hip and not outwards,
		hence the change in sign. This transform is that of the robot when it uses the reamer.
	*/

	auto vectorX = pImplantToBoneTransform->GetMatrix().GetVnlMatrix().get_column(0);
	Point implantVectorZ = -Point(vectorX[0], vectorX[1], vectorX[2]);
	implantVectorZ.normalice();

	Point refPoint = mPelvis.getCoronalPlaneAPP().getProjectionPoint(mPelvis.getPubicJoin());

	Point pelvisVectorAP = mPelvis.getPelvisVectorAP();

	/*
		The plane formed by the normal of the implant and the AP vector of the pelvis is used.
		The anteversion rotation occurs in that plane, therefore in that plane the angle is sought
	*/

	Point interceptionVector = mPelvis.getCoronalPlaneAPP().getInterceptionPlaneVector(implantVectorZ, pelvisVectorAP);
	interceptionVector.normalice();

	double angle1 = ImplantTools::getAngleBetweenVectorsDegree(interceptionVector, implantVectorZ);
	double angle2 = ImplantTools::getAngleBetweenVectorsDegree(-interceptionVector, implantVectorZ);
	double angle = std::min(angle1, angle2);

	refPoint = refPoint + 1000. * implantVectorZ;

	if (mPelvis.getCoronalPlaneAPP().eval(refPoint) < 0)
	{
		return angle;
	}
	else
	{
		return -angle;
	}
}

double HipPelvisImplantsMatchInfo::getCupShiftSuperior(const itk::Rigid3DTransform<>::Pointer pImplantToBoneTransform) const
{
	auto centerFromTransform = pImplantToBoneTransform->GetCenter();
	Point cupCenter = Point(centerFromTransform[0], centerFromTransform[1], centerFromTransform[2]);

	Line superiorLine(mPelvis.getPelvisVectorInfSup(), mHipCenterOfRotation);

	Plane axial;
	axial.init(mPelvis.getPelvisVectorInfSup(), mHipCenterOfRotation);

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

double HipPelvisImplantsMatchInfo::getCupShiftLateral(const itk::Rigid3DTransform<>::Pointer pImplantToBoneTransform) const
{
	auto centerFromTransform = pImplantToBoneTransform->GetCenter();
	Point cupCenter = Point(centerFromTransform[0], centerFromTransform[1], centerFromTransform[2]);

	Line lateralLine(mPelvis.getPelvisVectorLateralASIS(), mHipCenterOfRotation);

	Plane sagital;
	sagital.init(mPelvis.getPelvisVectorLateralASIS(), mHipCenterOfRotation);

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

double HipPelvisImplantsMatchInfo::getCupShiftAnterior(const itk::Rigid3DTransform<>::Pointer pImplantToBoneTransform) const
{
	auto centerFromTransform = pImplantToBoneTransform->GetCenter();
	Point cupCenter = Point(centerFromTransform[0], centerFromTransform[1], centerFromTransform[2]);

	Line anteriorLine(mPelvis.getPelvisVectorAP(), mHipCenterOfRotation);

	Plane coronal;
	coronal.init(mPelvis.getPelvisVectorAP(), mHipCenterOfRotation);

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

double HipPelvisImplantsMatchInfo::getCupInclination(const Point& pVectorToHipCenter) const
{
	Point rotationAxis = mPelvis.getPelvisVectorAP();
	Point referenceVector = mPelvis.getPelvisVectorInfSup();

	Point implantVector = pVectorToHipCenter;
	implantVector.normalice();

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

double HipPelvisImplantsMatchInfo::getCupAntversion(const Point& pVectorToHipCenter) const
{
	Point implantVector = pVectorToHipCenter;
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

double HipPelvisImplantsMatchInfo::getCupShiftSuperior(const Point& pCenterOfRotation) const
{
	Line superiorLine(mPelvis.getPelvisVectorInfSup(), mHipCenterOfRotation);

	Plane axial;
	axial.init(mPelvis.getPelvisVectorInfSup(), mHipCenterOfRotation);

	Point cupCenter = pCenterOfRotation;
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

double HipPelvisImplantsMatchInfo::getCupShiftLateral(const Point& pCenterOfRotation) const
{
	Line lateralLine(mPelvis.getPelvisVectorLateralASIS(), mHipCenterOfRotation);

	Plane sagital;
	sagital.init(mPelvis.getPelvisVectorLateralASIS(), mHipCenterOfRotation);

	Point cupCenter = pCenterOfRotation;
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

double HipPelvisImplantsMatchInfo::getCupShiftAnterior(const Point& pCenterOfRotation) const
{
	Line anteriorLine(mPelvis.getPelvisVectorAP(), mHipCenterOfRotation);

	Plane coronal;
	coronal.init(mPelvis.getPelvisVectorAP(), mHipCenterOfRotation);

	Point cupCenter = pCenterOfRotation;
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

double HipPelvisImplantsMatchInfo::getReamingDistance(const Point& pCenterOfRotation) const
{
	cv::Mat cupCenter = mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint() + mTranslationCup;
	return ImplantTools::getDistanceBetweenPoints(pCenterOfRotation, Point(cupCenter));
}

double HipPelvisImplantsMatchInfo::getStemAxisShiftSuperiorHip() const
{
	//cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	//stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	//Point stemCenter = Point(stemHeadCenter);

	cv::Mat stemBasePoint = (mRotationStem * mImplantStem.getBasePoint().ToMatPoint()) + mTranslationStem;
	Point stemBase = Point(stemBasePoint);

	Line superiorLine(mPelvis.getPelvisVectorInfSup(), mHipCenterOfRotation);

	Plane axial;
	axial.init(mPelvis.getPelvisVectorInfSup(), mHipCenterOfRotation);
	
	stemBase = superiorLine.getProjectPoint(stemBase);
	double distance = ImplantTools::getDistanceBetweenPoints(stemBase, mHipCenterOfRotation);

	if (axial.eval(stemBase) >= 0)
	{
		return distance;
	}
	else
	{
		return -distance;
	}
}

double HipPelvisImplantsMatchInfo::getStemAxisShiftLateralHip() const
{
	//cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	//stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	//Point stemCenter = Point(stemHeadCenter);

	cv::Mat stemBasePoint = (mRotationStem * mImplantStem.getBasePoint().ToMatPoint()) + mTranslationStem;
	Point stemBase = Point(stemBasePoint);

	Line lateralLine(mPelvis.getPelvisVectorLateralASIS(), mHipCenterOfRotation);

	Plane sagital;
	sagital.init(mPelvis.getPelvisVectorLateralASIS(), mHipCenterOfRotation);

	stemBase = lateralLine.getProjectPoint(stemBase);
	double distance = ImplantTools::getDistanceBetweenPoints(stemBase, mHipCenterOfRotation);

	if (sagital.eval(stemBase) >= 0)
	{
		return distance;
	}
	else
	{
		return -distance;
	}
}

double HipPelvisImplantsMatchInfo::getStemAxisShiftAnteriorHip() const
{
	//cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	//stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	//Point stemCenter = Point(stemHeadCenter);

	cv::Mat stemBasePoint = (mRotationStem * mImplantStem.getBasePoint().ToMatPoint()) + mTranslationStem;
	Point stemBase = Point(stemBasePoint);

	Line anteriorLine(mPelvis.getPelvisVectorAP(), mHipCenterOfRotation);

	Plane coronal;
	coronal.init(mPelvis.getPelvisVectorAP(), mHipCenterOfRotation);

	stemBase = anteriorLine.getProjectPoint(stemBase);
	double distance = ImplantTools::getDistanceBetweenPoints(stemBase, mHipCenterOfRotation);

	if (coronal.eval(stemBase) >= 0)
	{
		return distance;
	}
	else
	{
		return -distance;
	}
}

double HipPelvisImplantsMatchInfo::getStemAxisShiftSuperiorCup() const
{
	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;
	Point cupCenter = Point(cupCenterMat);

	//cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	//stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	//Point stemBase = Point(stemHeadCenter);

	cv::Mat stemBasePoint = (mRotationStem * mImplantStem.getBasePoint().ToMatPoint()) + mTranslationStem;
	Point stemBase = Point(stemBasePoint);

	Line superiorLine(mPelvis.getPelvisVectorInfSup(), cupCenter);

	Plane axial;
	axial.init(mPelvis.getPelvisVectorInfSup(), cupCenter);

	stemBase = superiorLine.getProjectPoint(stemBase);
	double distance = ImplantTools::getDistanceBetweenPoints(stemBase, cupCenter);

	if (axial.eval(stemBase) >= 0)
	{
		return distance;
	}
	else
	{
		return -distance;
	}
}

double HipPelvisImplantsMatchInfo::getStemAxisShiftLateralCup() const
{
	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;
	Point cupCenter = Point(cupCenterMat);

	//cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	//stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	//Point stemBase = Point(stemHeadCenter);

	cv::Mat stemBasePoint = (mRotationStem * mImplantStem.getBasePoint().ToMatPoint()) + mTranslationStem;
	Point stemBase = Point(stemBasePoint);

	Line lateralLine(mPelvis.getPelvisVectorLateralASIS(), cupCenter);

	Plane sagital;
	sagital.init(mPelvis.getPelvisVectorLateralASIS(), cupCenter);

	stemBase = lateralLine.getProjectPoint(stemBase);
	double distance = ImplantTools::getDistanceBetweenPoints(stemBase, cupCenter);

	if (sagital.eval(stemBase) >= 0)
	{
		return distance;
	}
	else
	{
		return -distance;
	}
}

double HipPelvisImplantsMatchInfo::getStemAxisShiftAnteriorCup() const
{
	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;
	Point cupCenter = Point(cupCenterMat);

	//cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	//stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	//Point stemBase = Point(stemHeadCenter);

	cv::Mat stemBasePoint = (mRotationStem * mImplantStem.getBasePoint().ToMatPoint()) + mTranslationStem;
	Point stemBase = Point(stemBasePoint);

	Line anteriorLine(mPelvis.getPelvisVectorAP(), cupCenter);

	Plane coronal;
	coronal.init(mPelvis.getPelvisVectorAP(), cupCenter);

	stemBase = anteriorLine.getProjectPoint(stemBase);
	double distance = ImplantTools::getDistanceBetweenPoints(stemBase, cupCenter);

	if (coronal.eval(stemBase) >= 0)
	{
		return distance;
	}
	else
	{
		return -distance;
	}
}

double HipPelvisImplantsMatchInfo::getStemHeadShiftSuperiorCup() const
{
	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;
	Point cupCenter = Point(cupCenterMat);

	cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	Point stemBase = Point(stemHeadCenter);

	Line superiorLine(mPelvis.getPelvisVectorInfSup(), cupCenter);

	Plane axial;
	axial.init(mPelvis.getPelvisVectorInfSup(), cupCenter);

	stemBase = superiorLine.getProjectPoint(stemBase);
	double distance = ImplantTools::getDistanceBetweenPoints(stemBase, cupCenter);

	if (axial.eval(stemBase) >= 0)
	{
		return distance;
	}
	else
	{
		return -distance;
	}
}

double HipPelvisImplantsMatchInfo::getStemHeadShiftLateralCup() const
{
	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;
	Point cupCenter = Point(cupCenterMat);

	cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	Point stemBase = Point(stemHeadCenter);

	Line lateralLine(mPelvis.getPelvisVectorLateralASIS(), cupCenter);

	Plane sagital;
	sagital.init(mPelvis.getPelvisVectorLateralASIS(), cupCenter);

	stemBase = lateralLine.getProjectPoint(stemBase);
	double distance = ImplantTools::getDistanceBetweenPoints(stemBase, cupCenter);

	if (sagital.eval(stemBase) >= 0)
	{
		return distance;
	}
	else
	{
		return -distance;
	}
}

double HipPelvisImplantsMatchInfo::getStemHeadShiftAnteriorCup() const
{
	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;
	Point cupCenter = Point(cupCenterMat);

	cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	Point stemBase = Point(stemHeadCenter);

	Line anteriorLine(mPelvis.getPelvisVectorAP(), cupCenter);

	Plane coronal;
	coronal.init(mPelvis.getPelvisVectorAP(), cupCenter);

	stemBase = anteriorLine.getProjectPoint(stemBase);
	double distance = ImplantTools::getDistanceBetweenPoints(stemBase, cupCenter);

	if (coronal.eval(stemBase) >= 0)
	{
		return distance;
	}
	else
	{
		return -distance;
	}
}

void HipPelvisImplantsMatchInfo::setStemTransform(const itk::Rigid3DTransform<>::Pointer pImplantToBoneStemTransform)
{
	mRotationStem = Rigid3DTransformToCVRotation(pImplantToBoneStemTransform);
	mTranslationStem = Rigid3DTransformToCVTranslation(pImplantToBoneStemTransform);
}

void HipPelvisImplantsMatchInfo::setStemHeadTransform(const itk::Rigid3DTransform<>::Pointer pImplantHeadToStemTransform)
{
	mRotationStemHead = Rigid3DTransformToCVRotation(pImplantHeadToStemTransform);
	mTranslationStemHead = Rigid3DTransformToCVTranslation(pImplantHeadToStemTransform);
}

itk::Rigid3DTransform<>::Pointer HipPelvisImplantsMatchInfo::getITKStemTransform() const
{
	itk::Matrix< double, 3, 3 > rotation;
	itk::Vector< double, 3 > translate;

	rotation[0][0] = mRotationStem.at<double>(0, 0);
	rotation[0][1] = mRotationStem.at<double>(0, 1);
	rotation[0][2] = mRotationStem.at<double>(0, 2);

	rotation[1][0] = mRotationStem.at<double>(1, 0);
	rotation[1][1] = mRotationStem.at<double>(1, 1);
	rotation[1][2] = mRotationStem.at<double>(1, 2);

	rotation[2][0] = mRotationStem.at<double>(2, 0);
	rotation[2][1] = mRotationStem.at<double>(2, 1);
	rotation[2][2] = mRotationStem.at<double>(2, 2);

	translate[0] = mTranslationStem.at<double>(0, 0);
	translate[1] = mTranslationStem.at<double>(1, 0);
	translate[2] = mTranslationStem.at<double>(2, 0);

	itk::Rigid3DTransform<double>::Pointer transform = itk::VersorRigid3DTransform<double>::New();
	transform->SetMatrix(rotation);
	transform->SetOffset(translate);

	return transform;
}

double HipPelvisImplantsMatchInfo::getStemVersion() const
{
	cv::Mat neckAxisMat = mRotationStem * mImplantStem.getVectorNeckToHead().ToMatPoint();
	Point neckAxis = Point(neckAxisMat);
	double degree = mPelvis.getFemurVersionDegree(neckAxis);
	return degree;
}

double HipPelvisImplantsMatchInfo::getRealStemVersion(const Point& pVectorNeckToHead) const
{
	double degree = mPelvis.getFemurVersionDegree(pVectorNeckToHead);
	return degree;
}

double HipPelvisImplantsMatchInfo::getCombinedOffsetDistance() const
{
	/*
	cv::Mat stemBasePointMat = (mRotationStem * mImplantStem.getBasePoint().ToMatPoint()) + mTranslationStem;
	Point stemBasePoint = Point(stemBasePointMat);

	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;
	Point cupCenter = Point(cupCenterMat);

	Point axisLeg = mPelvis.getFemurOperationSide().getCanalAxisVectorInfSup();
	Line axisLine(axisLeg, mPelvis.getFemurOperationSide().getCanalAxisPoint());
	Point translation = stemBasePoint - axisLine.getProjectPoint(stemBasePoint);
	translation = translation + (cupCenter - mPelvis.getFemurOperationSide().getHeadCenter());
	*/

	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;
	Point cupCenter = Point(cupCenterMat);

	cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	Point stemCenter = Point(stemHeadCenter);

	Point translation = stemCenter - cupCenter;


	return mPelvis.getCombinedOffsetDistance(cupCenter, translation.ToMatPoint());
}

double HipPelvisImplantsMatchInfo::getRealCombinedOffsetDistance(const Point& pFinalCupCenter) const
{
	cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	Point stemCenter = Point(stemHeadCenter);

	Point translation = stemCenter - pFinalCupCenter;

	return mPelvis.getCombinedOffsetDistance(pFinalCupCenter, translation.ToMatPoint());
}

double HipPelvisImplantsMatchInfo::getHipLengthDistance() const
{
	/*cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	Point stemCenter = Point(stemHeadCenter);

	cv::Mat stemBasePointMat = (mRotationStem * mImplantStem.getBasePoint().ToMatPoint()) + mTranslationStem;
	Point stemBasePoint = Point(stemBasePointMat);

	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;
	Point cupCenter = Point(cupCenterMat);

	Point axisLeg = mPelvis.getFemurOperationSide().getCanalAxisVectorInfSup();
	Line axisLine(axisLeg, mPelvis.getFemurOperationSide().getCanalAxisPoint());

	Point translation = stemBasePoint - axisLine.getProjectPoint(stemBasePoint);
	translation = translation + (cupCenter - mPelvis.getFemurOperationSide().getHeadCenter());

	Point hipOnAxis = axisLine.getProjectPoint(mPelvis.getFemurOperationSide().getHeadCenter());
	Point kneeOnAxis = axisLine.getProjectPoint(mPelvis.getFemurOperationSide().getKneeCenter());
	Point stemOnAxis = axisLine.getProjectPoint(stemCenter);

	double distance = ImplantTools::getDistanceBetweenPoints(hipOnAxis, kneeOnAxis);
	double distanceTemp = ImplantTools::getDistanceBetweenPoints(stemOnAxis, kneeOnAxis);
	double diff = distanceTemp - distance;
	axisLeg.normalice();

	if (diff > 0)
	{
		translation = translation + diff * axisLeg;
	}
	else
	{
		translation = translation + diff * axisLeg;
	}*/
	
	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;
	Point cupCenter = Point(cupCenterMat);

	cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	Point stemCenter = Point(stemHeadCenter);

	Point translation = stemCenter - cupCenter;

	return mPelvis.getHipLengthDistance(cupCenter, translation.ToMatPoint());
}

double HipPelvisImplantsMatchInfo::getRealHipLengthDistance(const Point& pFinalCupCenter) const
{
	cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	Point stemCenter = Point(stemHeadCenter);

	Point translation = stemCenter - pFinalCupCenter;

	return mPelvis.getHipLengthDistance(pFinalCupCenter, translation.ToMatPoint());
}

double HipPelvisImplantsMatchInfo::getCoverageFraction() const
{
	double area = mImplantCup.getHemiSphereSurfaceArea();
	auto hemiSphere = mImplantCup.getHemiSphereCup();

	auto cup = ImplantTools::transformPolydata(hemiSphere, mRotationCup, mTranslationCup);

	double overlappingArea = ImplantTools::getOverlappingArea(cup, mPelvis.getPelvisVTK(), mPelvis.getImplicitPelvisDistance());

	if (area > 0)
	{
		return overlappingArea / area;
	}

	return 0;
}

itk::Rigid3DTransform<>::Pointer HipPelvisImplantsMatchInfo::setCupAngles(double pAbductionAngle, double pAnteversionAngle)
{
	std::vector<cv::Point3d> implantVectors;
	std::vector<cv::Point3d> pelvisVectors;

	Point implantX = mImplantCup.getVectorX();
	Point implantZ = mImplantCup.getVectorZ();
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

	cv::Mat CenterBefore = mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint() + mTranslationCup;

	cv::Mat inverse = (implantMatrix.t()).inv();
	mRotationCup = (pelvisMatrix.t()) * inverse;

	cv::Mat CenterAfter = mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint() + mTranslationCup;

	cv::Mat diff = CenterAfter - CenterBefore;
	
	mTranslationCup = mTranslationCup - diff;

	return getITKCupTransform();
}

itk::Vector< double, 3 > HipPelvisImplantsMatchInfo::setCupTranslation(double pShifSuperior, double pShifLateral, double pShiftAnterior)
{
	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;
	Point cupCenter = Point(cupCenterMat);

	Point newCenterCup = mHipCenterOfRotation + pShifSuperior * mPelvis.getPelvisVectorInfSup();
	newCenterCup = newCenterCup + pShifLateral * mPelvis.getPelvisVectorLateralASIS();
	newCenterCup = newCenterCup + pShiftAnterior * mPelvis.getPelvisVectorAP();

	Point diff = newCenterCup - cupCenter;
	mTranslationCup = mTranslationCup + diff.ToMatPoint();

	return ImplantTools::CVTranslationToITKVector(mTranslationCup);
}

itk::Vector< double, 3 > HipPelvisImplantsMatchInfo::setStemHeadTranslation(double pShifSuperior, double pShifLateral, double pShiftAnterior)
{
	cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	Point stemCenter = Point(stemHeadCenter);

	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;
	Point cupCenter = Point(cupCenterMat);

	Point newCenterStem = cupCenter + pShifSuperior * mPelvis.getPelvisVectorInfSup();
	newCenterStem = newCenterStem + pShifLateral * mPelvis.getPelvisVectorLateralASIS();
	newCenterStem = newCenterStem + pShiftAnterior * mPelvis.getPelvisVectorAP();

	Point diff = newCenterStem - stemCenter;
	mTranslationStem = mTranslationStem + diff.ToMatPoint();

	return ImplantTools::CVTranslationToITKVector(mTranslationStem);
}

itk::Vector< double, 3 > HipPelvisImplantsMatchInfo::matchStemToHipRotationCenter()
{
	cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	Point stemCenter = Point(stemHeadCenter);

	Point diff = mHipCenterOfRotation - stemCenter;
	mTranslationStem = mTranslationStem + diff.ToMatPoint();

	return ImplantTools::CVTranslationToITKVector(mTranslationStem);
}

itk::Vector< double, 3 > HipPelvisImplantsMatchInfo::matchStemToCupRotationCenter()
{
	cv::Mat cupCenterMat = (mRotationCup * mImplantCup.getCenterOfRotationImplant().ToMatPoint()) + mTranslationCup;
	Point cupCenter = Point(cupCenterMat);

	cv::Mat stemHeadCenter = (mRotationStemHead * mImplantStemHead.getCenterOfSphere().ToMatPoint()) + mTranslationStemHead;
	stemHeadCenter = (mRotationStem * stemHeadCenter) + mTranslationStem;
	Point stemCenter = Point(stemHeadCenter);

	Point diff = cupCenter - stemCenter;
	mTranslationStem = mTranslationStem + diff.ToMatPoint();

	return ImplantTools::CVTranslationToITKVector(mTranslationStem);
}

itk::Rigid3DTransform<>::Pointer HipPelvisImplantsMatchInfo::setStemVersionAngle(double pStemVersionAngleDegree)
{
	cv::Mat neckAxisMat = mRotationStem * mImplantStem.getVectorNeckToHeadPerpendicularToInfSup().ToMatPoint();
	Point neckAxis = Point(neckAxisMat);

	double angle = mPelvis.getFemurVersionRadian(neckAxis);
	double refAngle = (pStemVersionAngleDegree * PI) / 180.;

	Point canalAxis = mPelvis.getFemurOperationSide().getCanalAxisVectorInfSup();

	cv::Mat transformation;

	if (refAngle == angle)
	{
		return getITKStemTransform();
	}
	else if (refAngle > angle)
	{
		double temp = refAngle - angle;

		if (mPelvis.getSide() == PelvisSide::RIGHT_SIDE)
		{
			transformation = ImplantTools::getRotateMatrix(-canalAxis, temp);
		}
		else
		{
			transformation = ImplantTools::getRotateMatrix(canalAxis, temp);
		}
	}
	else
	{
		double temp = angle - refAngle;

		if (mPelvis.getSide() == PelvisSide::RIGHT_SIDE)
		{
			transformation = ImplantTools::getRotateMatrix(canalAxis, temp);
		}
		else
		{
			transformation = ImplantTools::getRotateMatrix(-canalAxis, temp);
		}
	}

	Plane planeProj, planeSource;
	planeProj.init(canalAxis, mPelvis.getFemurOperationSide().getCanalAxisPoint());
	planeSource.init(mImplantStem.getVectorInfSup(), mImplantStem.getCanalAxisTopPoint());
	planeSource.transformPlane(mRotationStem, mTranslationStem);
	auto mainVectorMat = transformation * (planeProj.getProjectionVector(neckAxis).ToMatPoint());
	Point mainVector = Point(mainVectorMat);

	Point newVectorFromImplant = ImplantTools::getOriginalVectorFromProjectionWithPlanes(planeProj, mainVector, planeSource);

	double myAngle = ImplantTools::getAngleBetweenVectors(newVectorFromImplant, neckAxis);

	auto newRotation = ImplantTools::getRotateMatrix(neckAxis.cross(newVectorFromImplant), myAngle);

	cv::Mat stemBaseBefore = mRotationStem * mImplantStem.getBasePoint().ToMatPoint() + mTranslationStem;

	mRotationStem = newRotation * mRotationStem;

	cv::Mat stemBaseAfter = mRotationStem * mImplantStem.getBasePoint().ToMatPoint() + mTranslationStem;
	cv::Mat diff = stemBaseBefore - stemBaseAfter;
	mTranslationStem += diff;

	return getITKStemTransform();
}


Plane HipPelvisImplantsMatchInfo::TransformPlane(const Plane& plane, const cv::Mat& rotation, const cv::Mat& translation) const
{
	cv::Mat transformNormalVector = rotation * plane.getNormalVectorMat();
	cv::Mat transformPoint = (rotation * plane.getPointMat()) + translation;
	Plane transformPlane;
	transformPlane.init(Point(transformNormalVector), Point(transformPoint));
	return transformPlane;
}