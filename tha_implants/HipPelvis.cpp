#include "vtkPolyDataConnectivityFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkClipClosedSurface.h"
#include "vtkPlane.h"
#include "vtkPlaneCollection.h"
#include "vtkCutter.h"
#include "vtkAppendPolyData.h"
#include "HipPelvis.hpp"
#include "ImplantTools.hpp"
#include "ImplantsException.hpp"
#include <opencv2/opencv.hpp>
#include "vtkMassProperties.h"
#include "vtkSphereSource.h"

using namespace THA::IMPLANTS;

HipPelvis::HipPelvis()
{
    isInit = false;
}


void HipPelvis::init(const Point& pLeftASIS, const Point& pRightASIS, const Point& pLeftPubicTubercle, const Point& pRightPubicTubercle, const vtkSmartPointer<vtkPolyData>& pPelvis,
	const HipFemur& pFemur, const HipFemurOppside& pFemurOppside, PelvisSide pSide, const Point& pHipCenterOfRotation, const Plane& pCoronalCT)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_PELVIS;
    }

    /*vtkNew<vtkPolyDataConnectivityFilter> pelvisConnectivityFilter;
    pelvisConnectivityFilter->SetInputData(pPelvis);
    pelvisConnectivityFilter->SetExtractionModeToLargestRegion();
    pelvisConnectivityFilter->Update();

    vtkNew<vtkCleanPolyData> CleanPelvis;
    CleanPelvis->SetInputData(pelvisConnectivityFilter->GetOutput());
    CleanPelvis->Update();

    this->mPelvis = CleanPelvis->GetOutput();*/

    this->mPelvis = pPelvis;
    this->mLeftASIS = pLeftASIS;
    this->mRightASIS = pRightASIS;
    this->mLeftPubicTubercle = pLeftPubicTubercle;
    this->mRightPubicTubercle = pRightPubicTubercle;
    this->mPubicJoin = (pLeftPubicTubercle + pRightPubicTubercle) / 2.0;
    this->mSide = pSide;

    Point vector1 = mRightASIS - mLeftASIS;
    Point vector2 = mPubicJoin - mLeftASIS;
    Point normalRef = vector1.cross(vector2);

    mPlaneAPP.init(normalRef, mPubicJoin);

	mFemurOperationSide = pFemur;
	mFemurOppsite = pFemurOppside;

	mImplicitPelvisDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
	mImplicitPelvisDistance->SetInput(pPelvis);

	mHipCenterOfRotation = pHipCenterOfRotation;

	coronalTiltAngle = getCoronalTiltAngle(pCoronalCT);

    isInit = true;
}

HipPelvis HipPelvis::getHipPelvisCopyObj(double pTiltAngleDegree) const
{
	HipPelvis copyObj = *this;
	copyObj.setCoronalTiltAngleDegree(pTiltAngleDegree);
	return copyObj;
}

void HipPelvis::setHipCenterOfRotation(const Point& pHipCenterOfRotation)
{
	mHipCenterOfRotation = pHipCenterOfRotation;
}

Point HipPelvis::getHipCenterOfRotation() const
{
	return mHipCenterOfRotation;
}

//Point HipPelvis::getMidASIS() const
//{
//    Plane sagital;
//    sagital.init(getPelvisVectorASIS(), mPubicJoin);
//
//    Line lineASIS = Line::makeLineWithPoints(mRightASIS, mLeftASIS);
//    Point midASIS = sagital.getInterceptionLinePoint(lineASIS);
//    return midASIS;
//}

Point HipPelvis::getRightASIS() const
{
	return mRightASIS;
}

Point HipPelvis::getLeftASIS() const
{
	return mLeftASIS;
}

Point HipPelvis::getPelvisVectorASIS() const
{
    Point vector = mRightASIS - mLeftASIS;
    vector.normalice();
    return vector;
}

Point HipPelvis::getPelvisVectorAP() const
{
    return getCoronalPlaneAPP().getNormalVector();
}

Plane HipPelvis::getCoronalPlaneAPP() const
{
	if (coronalTiltAngle == 0)
	{
		return mPlaneAPP;
	}

	cv::Mat rotation = ImplantTools::getRotateMatrix(getPelvisVectorASIS(), -coronalTiltAngle);
	cv::Mat rotVector = rotation * mPlaneAPP.getNormalVectorMat();
	Plane tiltPlane;
	tiltPlane.init(Point(rotVector), mPubicJoin);
	return tiltPlane;
}

double HipPelvis::getCoronalTiltAngle(const Plane& pCoronalCT) const
{
	if (pCoronalCT.getIsInit() == false)
	{
		return 0;
	}

	Plane coronalCT, sagital;
	coronalCT.init(pCoronalCT.getNormalVector(), getPubicJoin());
	Point rightPoint = coronalCT.getProjectionPoint(mRightASIS);
	Point leftPoint = coronalCT.getProjectionPoint(mLeftASIS);

	Point vector1 = rightPoint - leftPoint;
	Point vector2 = mPubicJoin - leftPoint;
	Point normalRef = vector1.cross(vector2);

	coronalCT.reverseByNormal(normalRef);
	sagital.init(getPelvisVectorASIS(), getPubicJoin());
	Point coronalApVectorOnSagital = sagital.getProjectionVector(coronalCT.getNormalVector());
	coronalApVectorOnSagital.normalice();
	
	double angle = ImplantTools::getAngleBetweenVectors(mPlaneAPP.getNormalVector(), coronalApVectorOnSagital);

	if (angle == 0)
	{
		return angle;
	}

	Point crossVector = coronalApVectorOnSagital.cross(mPlaneAPP.getNormalVector());
	crossVector.normalice();

	if (crossVector.dot(getPelvisVectorASIS()) >= 0)
	{
		return -angle;
	}
	else
	{
		return angle;
	}
}

void HipPelvis::setCoronalTiltAngleDegree(double pTiltAngleDegree)
{
	coronalTiltAngle = pTiltAngleDegree * PI / 180;
}

double HipPelvis::getCoronalTiltAngleDegree() const
{
	return coronalTiltAngle * 180 / PI;
}

Point HipPelvis::getPelvisVectorInfSup() const
{
    Point vector = getPelvisVectorASIS().cross(getPelvisVectorAP());
    vector.normalice();
    return vector;
}

Point HipPelvis::getPelvisVectorLateralASIS() const
{
	if (mSide == PelvisSide::RIGHT_SIDE)
	{
		return getPelvisVectorASIS();
	}
	else
	{
		return -getPelvisVectorASIS();
	}
}

cv::Mat HipPelvis::getFemurMechanicalAlignmentRotation() const
{
	auto femurObj = mFemurOperationSide;

	Point axialPlaneNormal = getPelvisVectorInfSup();
	Point mechanicalAxis = femurObj.getHeadCenter() - femurObj.getKneeCenter();
	mechanicalAxis.normalice();
	auto rotation1 = ImplantTools::GetGeneralRotateTransformVectors(mechanicalAxis, axialPlaneNormal);

	/*Point refAP;
	if (mSide == PelvisSide::RIGHT_SIDE)
	{
		refAP = mechanicalAxis.cross(femurObj.getCanalAxisVectorInfSup());
	}
	else
	{
		refAP = femurObj.getCanalAxisVectorInfSup().cross(mechanicalAxis);
	}
	refAP.normalice();

	auto refAPMat = rotation1 * refAP.ToMatPoint();
	refAP = Point(refAPMat);
	refAP.normalice();

	auto rotation2 = ImplantTools::GetGeneralRotateTransformVectors(refAP, getPelvisVectorAP());
	auto rotation = rotation2 * rotation1;*/
	return rotation1;
}

cv::Mat HipPelvis::getFemurMechanicalAlignmentRotation(const Point& pCupCenter, const Point& pTranslation) const
{
	auto femurObj = mFemurOperationSide;

	Point axialPlaneNormal = getPelvisVectorInfSup();
	Point mechanicalAxis = pCupCenter - (femurObj.getKneeCenter() + pTranslation);
	mechanicalAxis.normalice();
	auto rotation1 = ImplantTools::GetGeneralRotateTransformVectors(mechanicalAxis, axialPlaneNormal);

	/*Point refAP;
	if (mSide == PelvisSide::RIGHT_SIDE)
	{
		refAP = mechanicalAxis.cross(femurObj.getCanalAxisVectorInfSup());
	}
	else
	{
		refAP = femurObj.getCanalAxisVectorInfSup().cross(mechanicalAxis);
	}
	refAP.normalice();

	auto refAPMat = rotation1 * refAP.ToMatPoint();
	refAP = Point(refAPMat);
	refAP.normalice();

	auto rotation2 = ImplantTools::GetGeneralRotateTransformVectors(refAP, getPelvisVectorAP());
	auto rotation = rotation2 * rotation1;*/
	return rotation1;
}

cv::Mat HipPelvis::getFemurMechanicalAlignmentRotationOppsite() const
{
	auto femurObj = mFemurOppsite;

	Point axialPlaneNormal = getPelvisVectorInfSup();
	Point mechanicalAxis = femurObj.getHeadCenter() - femurObj.getKneeCenter();
	mechanicalAxis.normalice();
	auto rotation1 = ImplantTools::GetGeneralRotateTransformVectors(mechanicalAxis, axialPlaneNormal);

	/*Point refAP;
	if (mSide == PelvisSide::RIGHT_SIDE)
	{
		refAP = femurObj.getCanalAxisVectorInfSup().cross(mechanicalAxis);
	}
	else
	{
		refAP = mechanicalAxis.cross(femurObj.getCanalAxisVectorInfSup());
	}
	refAP.normalice();

	auto refAPMat = rotation1 * refAP.ToMatPoint();
	refAP = Point(refAPMat);
	refAP.normalice();

	auto rotation2 = ImplantTools::GetGeneralRotateTransformVectors(refAP, getPelvisVectorAP());
	auto rotation = rotation2 * rotation1;*/
	return rotation1;
}

cv::Mat HipPelvis::getFemurCanalAlignmentRotation() const
{
	auto femurObj = mFemurOperationSide;

	Point axialPlaneNormal = getPelvisVectorInfSup();
	Point mechanicalAxis = femurObj.getHeadCenter() - femurObj.getKneeCenter();
	mechanicalAxis.normalice();
	auto rotation1 = ImplantTools::GetGeneralRotateTransformVectors(femurObj.getCanalAxisVectorInfSup(), axialPlaneNormal);

	/*Point refAP;
	if (mSide == PelvisSide::RIGHT_SIDE)
	{
		refAP = mechanicalAxis.cross(femurObj.getCanalAxisVectorInfSup());
	}
	else
	{
		refAP = femurObj.getCanalAxisVectorInfSup().cross(mechanicalAxis);
	}
	refAP.normalice();

	auto refAPMat = rotation1 * refAP.ToMatPoint();
	refAP = Point(refAPMat);
	refAP.normalice();

	auto rotation2 = ImplantTools::GetGeneralRotateTransformVectors(refAP, getPelvisVectorAP());
	auto rotation = rotation2 * rotation1;*/
	return rotation1;
}

cv::Mat HipPelvis::getFemurCanalAlignmentRotation(const Point& pCupCenter, const Point& pTranslation) const
{
	auto femurObj = mFemurOperationSide;

	Point axialPlaneNormal = getPelvisVectorInfSup();
	Point mechanicalAxis = pCupCenter - (femurObj.getKneeCenter() + pTranslation);
	mechanicalAxis.normalice();
	auto rotation1 = ImplantTools::GetGeneralRotateTransformVectors(femurObj.getCanalAxisVectorInfSup(), axialPlaneNormal);

	/*Point refAP;
	if (mSide == PelvisSide::RIGHT_SIDE)
	{
		refAP = mechanicalAxis.cross(femurObj.getCanalAxisVectorInfSup());
	}
	else
	{
		refAP = femurObj.getCanalAxisVectorInfSup().cross(mechanicalAxis);
	}
	refAP.normalice();

	auto refAPMat = rotation1 * refAP.ToMatPoint();
	refAP = Point(refAPMat);
	refAP.normalice();

	auto rotation2 = ImplantTools::GetGeneralRotateTransformVectors(refAP, getPelvisVectorAP());
	auto rotation = rotation2 * rotation1;*/
	return rotation1;
}

cv::Mat HipPelvis::getFemurCanalAlignmentRotationOppsite() const
{
	auto femurObj = mFemurOppsite;

	Point axialPlaneNormal = getPelvisVectorInfSup();
	Point mechanicalAxis = femurObj.getHeadCenter() - femurObj.getKneeCenter();
	mechanicalAxis.normalice();
	auto rotation1 = ImplantTools::GetGeneralRotateTransformVectors(femurObj.getCanalAxisVectorInfSup(), axialPlaneNormal);

	/*Point refAP;
	if (mSide == PelvisSide::RIGHT_SIDE)
	{
		refAP = femurObj.getCanalAxisVectorInfSup().cross(mechanicalAxis);
	}
	else
	{
		refAP = mechanicalAxis.cross(femurObj.getCanalAxisVectorInfSup());
	}
	refAP.normalice();

	auto refAPMat = rotation1 * refAP.ToMatPoint();
	refAP = Point(refAPMat);
	refAP.normalice();

	auto rotation2 = ImplantTools::GetGeneralRotateTransformVectors(refAP, getPelvisVectorAP());
	auto rotation = rotation2 * rotation1;*/
	return rotation1;
}

double HipPelvis::getHipLengthDistance() const
{
	auto femurObj = mFemurOperationSide;
	auto rotation = getFemurMechanicalAlignmentRotation();

	Point headCenter = femurObj.getHeadCenter();
	cv::Mat refMat = headCenter.ToMatPoint() + rotation * (femurObj.getLesserTrochanter().ToMatPoint() - headCenter.ToMatPoint());

	Point ref = mPlaneAPP.getProjectionPoint(Point(refMat));

	Line lineASIS = Line::makeLineWithPoints(mRightASIS, mLeftASIS);
	return lineASIS.getDistanceFromPoint(ref);
}

double HipPelvis::getHipLengthDistanceTest(const HipFemur& pFemurObj, const Plane& pPlaneAPP, const Point& PelvisVectorInfSup,
	const Point& pRightASIS, const Point& pLeftASIS) const
{
	auto femurObj = pFemurObj;
	Point axialPlaneNormal = PelvisVectorInfSup;
	Point mechanicalAxis = femurObj.getHeadCenter() - femurObj.getKneeCenter();
	mechanicalAxis.normalice();
	auto rotation1 = ImplantTools::GetGeneralRotateTransformVectors(mechanicalAxis, axialPlaneNormal);

	/*Point refAP;
	if (mSide == PelvisSide::RIGHT_SIDE)
	{
		refAP = mechanicalAxis.cross(femurObj.getCanalAxisVectorInfSup());
	}
	else
	{
		refAP = femurObj.getCanalAxisVectorInfSup().cross(mechanicalAxis);
	}
	refAP.normalice();

	auto refAPMat = rotation1 * refAP.ToMatPoint();
	refAP = Point(refAPMat);
	refAP.normalice();

	auto rotation2 = ImplantTools::GetGeneralRotateTransformVectors(refAP, getPelvisVectorAP());
	auto rotation = rotation2 * rotation1;*/

	Point headCenter = femurObj.getHeadCenter();
	cv::Mat refMat = headCenter.ToMatPoint() + rotation1 * (femurObj.getLesserTrochanter().ToMatPoint() - headCenter.ToMatPoint());
	Point ref = pPlaneAPP.getProjectionPoint(Point(refMat));

	Line lineASIS = Line::makeLineWithPoints(pRightASIS, pLeftASIS);
	return lineASIS.getDistanceFromPoint(ref);
}

double HipPelvis::getHipLengthDistanceOppsite() const
{
	auto femurObj = mFemurOppsite;
	auto rotation = getFemurMechanicalAlignmentRotationOppsite();

	Point headCenter = femurObj.getHeadCenter();

	cv::Mat refMat = headCenter.ToMatPoint() + rotation * (femurObj.getLesserTrochanter().ToMatPoint() - headCenter.ToMatPoint());
	Point ref = mPlaneAPP.getProjectionPoint(refMat);

	Line lineASIS = Line::makeLineWithPoints(mRightASIS, mLeftASIS);
	return lineASIS.getDistanceFromPoint(ref);
}


double HipPelvis::getCombinedOffsetDistance() const
{
	auto femurObj = mFemurOperationSide;

	auto rotation = getFemurCanalAlignmentRotation();
	Point headCenter = femurObj.getHeadCenter();

	cv::Mat canalPointMat = headCenter.ToMatPoint() + rotation * (femurObj.getCanalAxisPoint().ToMatPoint() - headCenter.ToMatPoint());
	Point canalPoint = mPlaneAPP.getProjectionPoint(Point(canalPointMat));

	Line midLine(getPelvisVectorInfSup(), mPubicJoin);
	
	return midLine.getDistanceFromPoint(canalPoint);
}

double HipPelvis::getCombinedOffsetDistanceOppsite() const
{
	auto femurObj = mFemurOppsite;
	auto rotation = getFemurCanalAlignmentRotationOppsite();

	Point headCenter =femurObj.getHeadCenter();

	cv::Mat canalPointMat = headCenter.ToMatPoint() + rotation * (femurObj.getCanalAxisPoint().ToMatPoint() - headCenter.ToMatPoint());
	Point canalPoint = mPlaneAPP.getProjectionPoint(Point(canalPointMat));

	Line midLine(getPelvisVectorInfSup(), mPubicJoin);

	return midLine.getDistanceFromPoint(canalPoint);
}

double HipPelvis::getHipLengthDistance(const Point& pHipHeadCenter, const cv::Mat& pFemurTranslation) const
{
	auto femurObj = mFemurOperationSide;

	Point femurMove = Point(pFemurTranslation);

	/*Point axialPlaneNormal = getPelvisVectorInfSup();
	Point mechanicalAxis = pHipHeadCenter - (femurObj.getKneeCenter() + femurMove);
	mechanicalAxis = mPlaneAPP.getProjectionVector(mechanicalAxis);
	mechanicalAxis.normalice();
	auto rotation = ImplantTools::GetGeneralRotateTransformVectors(mechanicalAxis, axialPlaneNormal);*/

	auto rotation = getFemurMechanicalAlignmentRotation(pHipHeadCenter, femurMove);
	cv::Mat refMat = pHipHeadCenter.ToMatPoint() + rotation * ((femurObj.getLesserTrochanter().ToMatPoint() + pFemurTranslation) - pHipHeadCenter.ToMatPoint());
	Point ref = mPlaneAPP.getProjectionPoint(Point(refMat));

	////////////////////////////////////////////////////////////////////////
	
	/*Point head = pHipHeadCenter;
	Point neckCenter = femurObj.getNeckCenter() + femurMove;
	Point canalCenter = femurObj.getCanalAxisPoint() + femurMove;
	Point trochanter = femurObj.getLesserTrochanter() + femurMove;
	Point medialEpi = femurObj.getMedialEpicondyle() + femurMove;
	Point lateralEpi = femurObj.getLateralEpicondyle() + femurMove;
	Point kneeCenter = femurObj.getKneeCenter() + femurMove;


	HipFemur pFemurObj;
	pFemurObj.init(head, neckCenter, canalCenter, trochanter, medialEpi, lateralEpi, kneeCenter, femurObj.getFemur());
	double result = getHipLengthDistanceTest(pFemurObj, mPlaneAPP, getPelvisVectorInfSup(), mRightASIS, mLeftASIS);

	std::cout << "Hip Lenght after move: " << result << std::endl;*/
	
	//////////////////////////////////////////////////////////////////////

	Line lineASIS = Line::makeLineWithPoints(mRightASIS, mLeftASIS);
	return lineASIS.getDistanceFromPoint(ref);
}

double HipPelvis::getCombinedOffsetDistance(const Point& pHipHeadCenter, const cv::Mat& pFemurTranslation) const
{
	auto femurObj = mFemurOperationSide;
	/*Point canalAxisVector = mPlaneAPP.getProjectionVector(femurObj.getCanalAxisVectorInfSup());
	canalAxisVector.normalice();
	Point axialPlaneNormal = getPelvisVectorInfSup();

	auto rotation = ImplantTools::GetGeneralRotateTransformVectors(canalAxisVector, axialPlaneNormal);*/

	Point femurMove = Point(pFemurTranslation);

	auto rotation = getFemurCanalAlignmentRotation(pHipHeadCenter, femurMove);

	cv::Mat canalPointMat = pHipHeadCenter.ToMatPoint() + rotation * ((femurObj.getCanalAxisPoint().ToMatPoint() + pFemurTranslation) - pHipHeadCenter.ToMatPoint());
	Point canalPoint = mPlaneAPP.getProjectionPoint(Point(canalPointMat));

	//Point centerASIS = (getRightASIS() + getLeftASIS()) / 2.;
	Line midLine(getPelvisVectorInfSup(), mPubicJoin);

	return midLine.getDistanceFromPoint(canalPoint);
}

double HipPelvis::getHipLengthDistance(const Point& pHipHeadCenter) const
{
	auto femurObj = mFemurOperationSide;

	Point femurMove = pHipHeadCenter - femurObj.getHeadCenter(); // check this!!!!

	/*Point axialPlaneNormal = getPelvisVectorInfSup();
	Point mechanicalAxis = pHipHeadCenter - (femurObj.getKneeCenter() + femurMove);
	mechanicalAxis = mPlaneAPP.getProjectionVector(mechanicalAxis);
	mechanicalAxis.normalice();
	auto rotation = ImplantTools::GetGeneralRotateTransformVectors(mechanicalAxis, axialPlaneNormal);*/

	auto rotation = getFemurMechanicalAlignmentRotation(pHipHeadCenter, femurMove);

	cv::Mat refMat = pHipHeadCenter.ToMatPoint() + rotation * ((femurObj.getLesserTrochanter().ToMatPoint() + femurMove.ToMatPoint()) - pHipHeadCenter.ToMatPoint());
	Point ref = mPlaneAPP.getProjectionPoint(Point(refMat));

	Line lineASIS = Line::makeLineWithPoints(mRightASIS, mLeftASIS);
	return lineASIS.getDistanceFromPoint(ref);
}

double HipPelvis::getCombinedOffsetDistance(const Point& pHipHeadCenter) const
{
	auto femurObj = mFemurOperationSide;
	/*Point canalAxisVector = mPlaneAPP.getProjectionVector(femurObj.getCanalAxisVectorInfSup());
	canalAxisVector.normalice();
	Point axialPlaneNormal = getPelvisVectorInfSup();

	auto rotation = ImplantTools::GetGeneralRotateTransformVectors(canalAxisVector, axialPlaneNormal);*/

	Point femurMove = pHipHeadCenter - femurObj.getHeadCenter();

	auto rotation = getFemurCanalAlignmentRotation(pHipHeadCenter, femurMove);

	cv::Mat canalPointMat = pHipHeadCenter.ToMatPoint() + rotation * ((femurObj.getCanalAxisPoint().ToMatPoint() + femurMove.ToMatPoint()) - pHipHeadCenter.ToMatPoint());
	Point canalPoint = mPlaneAPP.getProjectionPoint(Point(canalPointMat));

	//Point centerASIS = (getRightASIS() + getLeftASIS()) / 2.;
	Line midLine(getPelvisVectorInfSup(), mPubicJoin);

	return midLine.getDistanceFromPoint(canalPoint);
}

/*
double HipPelvis::getHipLengthDistance(const cv::Mat& pRotation, const cv::Mat& pTranslation) const
{
	auto femurObj = mFemurOperationSide;
	cv::Mat trochanterMat = pRotation * femurObj.getLesserTrochanter().ToMatPoint() + pTranslation;
	Point ref = mPlaneAPP.getProjectionPoint(Point(trochanterMat));

	Line lineASIS = Line::makeLineWithPoints(mRightASIS, mLeftASIS);
	return lineASIS.getDistanceFromPoint(ref);
}

double HipPelvis::getCombinedOffsetDistance(const cv::Mat& pRotation, const cv::Mat& pTranslation) const
{
	auto femurObj = mFemurOperationSide;

	cv::Mat canalAxisMat = pRotation * femurObj.getCanalAxisVectorInfSup().ToMatPoint();
	cv::Mat canalPointMat = pRotation * femurObj.getCanalAxisPoint().ToMatPoint() + pTranslation;
	cv::Mat kneeCenterMat = pRotation * femurObj.getKneeCenter().ToMatPoint() + pTranslation;
	Point kneeCenter = Point(kneeCenterMat);
	
	Point canalAxis = mPlaneAPP.getProjectionVector(canalAxisMat);
	Point canalPoint = mPlaneAPP.getProjectionPoint(canalPointMat);
	Line canalLine(canalAxis, canalPoint);
	Point kneeCenterOnCanalAxis = canalLine.getProjectPoint(kneeCenter);

	Point centerASIS = (getRightASIS() + getLeftASIS()) / 2.;
	Line midLine(getPelvisVectorInfSup(), centerASIS);

	return midLine.getDistanceFromPoint(kneeCenterOnCanalAxis);
}
*/

/*
std::pair<Point, Point> HipPelvis::getAbductionAnteversionVectorsZX(const Point& pCenterOfRotation, double pAbductionAngle, double pAnteversionAngle) const
{
	
	//The vector that goes from the base of the cup to its highest point is the one used as a reference to make the rotations.

	////////////////////////////////////////////// Abduction or incination

    Plane sagital;
    sagital.init(getPelvisVectorASIS(), mPubicJoin);
    sagital.reverseByPoint(pCenterOfRotation, false);

    Point rotateVector = getPelvisVectorInfSup();
    Point referenceVector = sagital.getNormalVector(); //rotateVector will rotate towards referenceVector.
    Point rotationAxis = rotateVector.cross(referenceVector); //it needs make cross because rotateVector will rotate towards referenceVector. RotationAxis is AP

    double angle = (pAbductionAngle * PI) / 180.0;
    cv::Mat rotMatrix = ImplantTools::getRotateMatrix(rotationAxis, angle);

    cv::Mat resultMat = rotMatrix * rotateVector.ToMatPoint();
    Point resultAbduction = Point(resultMat);
    resultAbduction.normalice();

    ////////////////////////////////////////////// Anteversion

    rotateVector = getPelvisVectorAP();
	referenceVector = resultAbduction;
	rotationAxis = rotateVector.cross(referenceVector);
	
	angle = (pAnteversionAngle * PI) / 180.0;
    rotMatrix = ImplantTools::getRotateMatrix(rotationAxis, angle);

    resultMat = rotMatrix * rotateVector.ToMatPoint();
    Point resultPointAnteversion = Point(resultMat);
    resultPointAnteversion.normalice();

	/////////////////////////////////// Last Rotation
	resultMat = rotMatrix * resultAbduction.ToMatPoint();
	resultAbduction = Point(resultMat);
	resultAbduction.normalice();

    return std::make_pair(resultAbduction, resultPointAnteversion);
}
*/

Point HipPelvis::getAbductionAnteversionVectorRotate(const Plane& pSagital, const Plane& pCoronal, const Point& pCenterOfRotation, double pAbductionAngle, double pAnteversionAngle) const
{
	/*
		The vector that goes from the base of the cup to its highest point is the one used as a reference to make the rotations.
	*/

	////////////////////////////////////////////// Abduction or inclination

	Plane sagital, coronal, axial;
	sagital = pSagital.getCopy();
	coronal = pCoronal.getCopy();
	axial.init(sagital.getNormalVector().cross(coronal.getNormalVector()), mFemurOperationSide.getKneeCenter());

	axial.reverseByPoint(mPubicJoin);
	sagital.movePlane(mRightASIS);
	sagital.reverseByPoint(mLeftASIS);
	coronal.reverseByNormal(axial.getNormalVector().cross(sagital.getNormalVector()));
	coronal.movePlane(getPubicJoin());

	Point rotateVector = axial.getNormalVector();
	Point rotationAxis = coronal.getNormalVector(); //RotationAxis is AP

	double angle = (pAbductionAngle * PI) / 180.0;
	cv::Mat rotMatrix = ImplantTools::getRotateMatrix(rotationAxis, angle);

	cv::Mat resultMat = rotMatrix * rotateVector.ToMatPoint();
	Point resultAbduction = Point(resultMat);
	resultAbduction.normalice();

	////////////////////////////////////////////// Anteversion
	rotateVector = resultAbduction;
	rotationAxis = coronal.getNormalVector().cross(rotateVector);
	angle = (pAnteversionAngle * PI) / 180.0;
	rotMatrix = ImplantTools::getRotateMatrix(rotationAxis, angle);
	double controlAngle = angle;
	double epsilon = 1e-9;
	for (int i = 0; i < 10; i++)
	{
		cv::Mat moveVector = rotMatrix * rotateVector.ToMatPoint();
		auto projVector = sagital.getProjectionVector(Point(moveVector));
		double tempAngle = ImplantTools::getAngleBetweenVectors(projVector, axial.getNormalVector());

		Point ref = getPubicJoin() + 1000. * projVector;

		if (coronal.eval(ref) < 0)
		{
			tempAngle = -tempAngle;
		}

		double diff = angle - tempAngle;
		if (std::abs(diff) > epsilon)
		{
			controlAngle += controlAngle * diff / tempAngle;
			rotMatrix = ImplantTools::getRotateMatrix(rotationAxis, controlAngle);
		}
		else
		{
			break;
		}
	}

	resultMat = rotMatrix * resultAbduction.ToMatPoint();
	resultAbduction = Point(resultMat);
	resultAbduction.normalice();
	return resultAbduction;
}

std::pair<cv::Point3d, double> HipPelvis::getNativeCenterOfRotation(const std::vector<Point>& pPoints)
{
    std::pair<cv::Point3d, double> result = ImplantTools::fitSphere(pPoints);
    if (result.second == -1)
    {
        throw ImplantExceptionCode::FAILED_COMPUTING_CENTER_OF_ROTATION_ON_HIP_PELVIS;
    }
    return result;
}

Point HipPelvis::getPubicJoin() const
{
    return mPubicJoin;
}

vtkSmartPointer<vtkPolyData> HipPelvis::getPelvisVTK() const
{
    return mPelvis;
}

PelvisSide HipPelvis::getSide() const
{
	return mSide;
}

double HipPelvis::getFemurVersionRadian(const Point& pNeckAxisVectorToHead) const
{
	return mFemurOperationSide.getFemurVersion(pNeckAxisVectorToHead, mSide);
}

double HipPelvis::getFemurVersionDegree(const Point& pNeckAxisVectorToHead) const
{
	double rad = mFemurOperationSide.getFemurVersion(pNeckAxisVectorToHead, mSide);
	return (rad * 180.) / PI;
}

double HipPelvis::getFemurVersionRadian() const
{
	return mFemurOperationSide.getFemurVersion(mFemurOperationSide.getNeckAxisVectorToHead(), mSide);
}

double HipPelvis::getFemurVersionDegree() const
{
	double rad = mFemurOperationSide.getFemurVersion(mFemurOperationSide.getNeckAxisVectorToHead(), mSide);
	return (rad * 180.) / PI;
}

HipFemur HipPelvis::getFemurOperationSide() const
{
	return mFemurOperationSide;
}

HipFemurOppside HipPelvis::getFemurOppsite() const
{
	return mFemurOppsite;
}

vtkSmartPointer<vtkImplicitPolyDataDistance> HipPelvis::getImplicitPelvisDistance() const
{
	return mImplicitPelvisDistance;
}

double HipPelvis::getNeckShaftAngle() const
{
	Point neckVector = mFemurOperationSide.getNeckAxisVectorToHead();
	Point canalVector = -mFemurOperationSide.getCanalAxisVectorInfSup();
	
	neckVector = mPlaneAPP.getProjectionVector(neckVector);
	canalVector = mPlaneAPP.getProjectionVector(canalVector);

	return ImplantTools::getAngleBetweenVectorsDegree(neckVector, canalVector);
}


/*
void HipPelvis::extractFemurAxisVector(bool surgerySide)
{
    Point vectorInfoSup = getPelvisVectorInfSup();

    Plane axial, sagital;
	Point initAxialPoint;
    Point refLesserTrochanter;

	if (surgerySide == true)
	{
		if (mSide == LEFT_SIDE)
		{
			refLesserTrochanter = mLeftLesserTrochanter;
			initAxialPoint = mLeftLesserTrochanter - 20 * vectorInfoSup;
		}
		else
		{
			refLesserTrochanter = mRightLesserTrochanter;
			initAxialPoint = mRightLesserTrochanter - 20 * vectorInfoSup;
		}
	}
	else
	{
		if (mSide == LEFT_SIDE)
		{
			refLesserTrochanter = mRightLesserTrochanter;
			initAxialPoint = mRightLesserTrochanter - 20 * vectorInfoSup;
		}
		else
		{
			refLesserTrochanter = mLeftLesserTrochanter;
			initAxialPoint = mLeftLesserTrochanter - 20 * vectorInfoSup;
		}
	}

    axial.init(vectorInfoSup, initAxialPoint);
    sagital.init(getPelvisVectorASIS(), mPubicJoin);

    axial.reverseByPoint(mLeftASIS, false);
    sagital.reverseByPoint(refLesserTrochanter);

    vtkNew<vtkPlane> vtkPlaneAxial, vtkPlaneSagital;

    Point myPoint = axial.getPoint();
    Point myNormal = axial.getNormalVector();
    vtkPlaneAxial->SetOrigin(myPoint.x, myPoint.y, myPoint.z);
    vtkPlaneAxial->SetNormal(myNormal.x, myNormal.y, myNormal.z);

    myPoint = sagital.getPoint();
    myNormal = sagital.getNormalVector();
    vtkPlaneSagital->SetOrigin(myPoint.x, myPoint.y, myPoint.z);
    vtkPlaneSagital->SetNormal(myNormal.x, myNormal.y, myNormal.z);

    vtkNew<vtkPlaneCollection> allPlanes;
    allPlanes->AddItem(vtkPlaneAxial);
    allPlanes->AddItem(vtkPlaneSagital);

    vtkNew<vtkClipClosedSurface> Clipper;
    Clipper->SetInputData(mPelvis);
    Clipper->SetClippingPlanes(allPlanes);
    Clipper->Update();

    auto femurData = Clipper->GetOutput();

    Point vectorNormalCut = vectorInfoSup;
    Point endAxialPoint;
    bool isOk = false;
    for (float i = 105; i >= 5; i -= 10)
    {
        Point myPointCut = initAxialPoint - i * vectorNormalCut;
        auto contour = ImplantTools::getContours(femurData, vectorNormalCut, myPointCut);
        if (contour->GetNumberOfPoints() != 0)
        {
            endAxialPoint = myPointCut;
            isOk = true;
            break;
        }
    }

    if (isOk == false)
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_FEMUR_AXIS_ON_HIP_PELVIS;
    }

    double femurSize = ImplantTools::getDistanceBetweenPoints(initAxialPoint, endAxialPoint);
    Point midFemurPoint = (initAxialPoint + endAxialPoint) / 2.0;
    std::vector<cv::Point3d> linearPoints;

    Point tAverage;
    for (float i = -0.2; i < 0.31; i += 0.1)
    {
        Point ringPoint = midFemurPoint + i * femurSize * vectorInfoSup;
        auto ringContour = ImplantTools::getMaxContour(femurData, vectorInfoSup, ringPoint);

        if (ringContour->GetNumberOfPoints() == 0)
        {
            throw ImplantExceptionCode::CAN_NOT_DETERMINE_FEMUR_AXIS_ON_HIP_PELVIS;
        }

        Point tCenter;
        ImplantTools::fitEllipse(ringContour, vectorInfoSup, tCenter);
        linearPoints.push_back(tCenter);

        tAverage = tAverage + tCenter;
    }

    tAverage = tAverage / double(linearPoints.size());

    cv::Vec6f outPut;
    cv::fitLine(linearPoints, outPut, CV_DIST_L2, 0, 0.01, 0.01);
    Point myVector = Point(outPut[0], outPut[1], outPut[2]);

    Line myLine(myVector, Point(outPut[3], outPut[4], outPut[5]));
    
    Plane planeTemp;
    planeTemp.init(myVector, endAxialPoint);
    planeTemp.reverseByPoint(initAxialPoint);

	if (surgerySide == true)
	{
		mFemurPointInsideCenter = myLine.getProjectPoint(tAverage);
		mFemurAxisVector = planeTemp.getNormalVector();
	}
	else
	{
		mFemurPointInsideCenterOppsite = myLine.getProjectPoint(tAverage);
		mFemurAxisVectorOppsite = planeTemp.getNormalVector();
	}
}
*/

vtkSmartPointer<vtkPolyData> HipPelvis::getPelvisFemurCutTest()
{
    Point vectorInfoSup = getPelvisVectorInfSup();

    Plane axial, sagital;
    Point initAxialPoint = ((mFemurOperationSide.getLesserTrochanter() + mFemurOppsite.getLesserTrochanter()) / 2.0) - 20 * vectorInfoSup;
	Point refLesserTrochanter = mFemurOperationSide.getLesserTrochanter();

    axial.init(vectorInfoSup, initAxialPoint);
    sagital.init(getPelvisVectorASIS(), mPubicJoin);

    axial.reverseByPoint(mLeftASIS, false);
    sagital.reverseByPoint(refLesserTrochanter);

    vtkNew<vtkPlane> vtkPlaneAxial, vtkPlaneSagital;

    Point myPoint = axial.getPoint();
    Point myNormal = axial.getNormalVector();
    vtkPlaneAxial->SetOrigin(myPoint.x, myPoint.y, myPoint.z);
    vtkPlaneAxial->SetNormal(myNormal.x, myNormal.y, myNormal.z);

    myPoint = sagital.getPoint();
    myNormal = sagital.getNormalVector();
    vtkPlaneSagital->SetOrigin(myPoint.x, myPoint.y, myPoint.z);
    vtkPlaneSagital->SetNormal(myNormal.x, myNormal.y, myNormal.z);

    vtkNew<vtkPlaneCollection> allPlanes;
    allPlanes->AddItem(vtkPlaneAxial);
    allPlanes->AddItem(vtkPlaneSagital);

    vtkNew<vtkClipClosedSurface> Clipper;
    Clipper->SetInputData(mPelvis);
    Clipper->SetClippingPlanes(allPlanes);
    Clipper->Update();

    auto femurData = Clipper->GetOutput();
    return femurData;
}