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
	const HipFemur& pFemur, const HipFemurOppside& pFemurOppside, PelvisSide pSide)
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

    isInit = true;
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
    return mPlaneAPP.getNormalVector();
}

Plane HipPelvis::getCoronalPlaneAPP() const
{
	return mPlaneAPP;
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

double HipPelvis::getHipLengthDistance() const
{
	auto femurObj = mFemurOperationSide;

	Point mechanicalAxis = femurObj.getHeadCenter() - femurObj.getKneeCenter();
	mechanicalAxis.normalice();
	mechanicalAxis = mPlaneAPP.getProjectionVector(mechanicalAxis);
	mechanicalAxis.normalice();
	Point refVector = getPelvisVectorInfSup();
	cv::Mat rotation = ImplantTools::GetGeneralRotateTransformVectors(mechanicalAxis, refVector);

	Point ref = mPlaneAPP.getProjectionPoint(femurObj.getLesserTrochanter());
	auto refMat = rotation * ref.ToMatPoint();
	ref = Point(refMat);

	Line lineASIS = Line::makeLineWithPoints(mRightASIS, mLeftASIS);
	return lineASIS.getDistanceFromPoint(ref);
}

double HipPelvis::getHipLengthDistanceOppsite() const
{
	auto femurObj = mFemurOppsite;

	Point mechanicalAxis = femurObj.getHeadCenter() - femurObj.getKneeCenter();
	mechanicalAxis.normalice();
	mechanicalAxis = mPlaneAPP.getProjectionVector(mechanicalAxis);
	mechanicalAxis.normalice();
	Point refVector = getPelvisVectorInfSup();
	cv::Mat rotation = ImplantTools::GetGeneralRotateTransformVectors(mechanicalAxis, refVector);

	Point ref = mPlaneAPP.getProjectionPoint(femurObj.getLesserTrochanter());
	auto refMat = rotation * ref.ToMatPoint();
	ref = Point(refMat);

	Line lineASIS = Line::makeLineWithPoints(mRightASIS, mLeftASIS);
	return lineASIS.getDistanceFromPoint(ref);
}

double HipPelvis::getCombinedOffsetDistance() const
{
	auto femurObj = mFemurOperationSide;
	Point canalAxis = mPlaneAPP.getProjectionVector(femurObj.getCanalAxisVectorInfSup());
	Point canalPoint = mPlaneAPP.getProjectionPoint(femurObj.getCanalAxisPoint());
	Line canalLine(canalAxis, canalPoint);
	Point kneeCenterOnCanalAxis = canalLine.getProjectPoint(femurObj.getKneeCenter());

	Point centerASIS = (getRightASIS() + getLeftASIS()) / 2.;
	Line midLine(getPelvisVectorInfSup(), centerASIS);
	
	return midLine.getDistanceFromPoint(kneeCenterOnCanalAxis);
}

double HipPelvis::getCombinedOffsetDistanceOppsite() const
{
	auto femurObj = mFemurOppsite;
	Point canalAxis = mPlaneAPP.getProjectionVector(femurObj.getCanalAxisVectorInfSup());
	Point canalPoint = mPlaneAPP.getProjectionPoint(femurObj.getCanalAxisPoint());
	Line canalLine(canalAxis, canalPoint);
	Point kneeCenterOnCanalAxis = canalLine.getProjectPoint(femurObj.getKneeCenter());

	Point centerASIS = (getRightASIS() + getLeftASIS()) / 2.;
	Line midLine(getPelvisVectorInfSup(), centerASIS);

	return midLine.getDistanceFromPoint(kneeCenterOnCanalAxis);
}

double HipPelvis::getHipLengthDistance(const Point& pMechanicalAxis) const
{
	auto femurObj = mFemurOperationSide;

	Point mechanicalAxis = pMechanicalAxis;
	mechanicalAxis.normalice();
	mechanicalAxis = mPlaneAPP.getProjectionVector(mechanicalAxis);
	mechanicalAxis.normalice();
	Point refVector = getPelvisVectorInfSup();
	cv::Mat rotation = ImplantTools::GetGeneralRotateTransformVectors(mechanicalAxis, refVector);

	Point ref = mPlaneAPP.getProjectionPoint(femurObj.getLesserTrochanter());
	auto refMat = rotation * ref.ToMatPoint();
	ref = Point(refMat);

	Line lineASIS = Line::makeLineWithPoints(mRightASIS, mLeftASIS);
	return lineASIS.getDistanceFromPoint(ref);
}

double HipPelvis::getCombinedOffsetDistance(const Point& pCanalAxis, const Point& pCanalAxisPoint) const
{
	auto femurObj = mFemurOperationSide;
	
	Point canalAxis = mPlaneAPP.getProjectionVector(pCanalAxis);
	Point canalPoint = mPlaneAPP.getProjectionPoint(pCanalAxisPoint);
	Line canalLine(canalAxis, canalPoint);
	Point kneeCenterOnCanalAxis = canalLine.getProjectPoint(femurObj.getKneeCenter());

	Point centerASIS = (getRightASIS() + getLeftASIS()) / 2.;
	Line midLine(getPelvisVectorInfSup(), centerASIS);

	return midLine.getDistanceFromPoint(kneeCenterOnCanalAxis);
}

std::pair<Point, Point> HipPelvis::getAbductionAnteversionVectorsZX(const Point& pCenterOfRotation, double pAbductionAngle, double pAnteversionAngle) const
{
	/*
		The vector that goes from the base of the cup to its highest point is the one used as a reference to make the rotations.
	*/

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

Point HipPelvis::getNativeCenterOfRotation(const std::vector<Point>& pPoints)
{
    std::pair<cv::Point3d, double> result = ImplantTools::fitSphere(pPoints);
    if (result.second == -1)
    {
        throw ImplantExceptionCode::FAILED_COMPUTING_CENTER_OF_ROTATION_ON_HIP_PELVIS;
    }
    return result.first;
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