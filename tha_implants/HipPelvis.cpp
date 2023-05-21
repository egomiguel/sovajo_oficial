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

#include "vtkSphereSource.h"

using namespace THA::IMPLANTS;

HipPelvis::HipPelvis()
{
    isInit = false;
}


void HipPelvis::init(const Point& pLeftASIS, const Point& pRightASIS, const Point& pLeftPubicTubercle, const Point& pRightPubicTubercle, const Point& pLeftLesserTrochanter, const Point& pRightLesserTrochanter, const vtkSmartPointer<vtkPolyData>& pPelvis, PelvisSide pSide)
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
    this->mLeftLesserTrochanter = pLeftLesserTrochanter;
    this->mRightLesserTrochanter = pRightLesserTrochanter;
    this->mPubicJoin = (pLeftPubicTubercle + pRightPubicTubercle) / 2.0;
    this->mSide = pSide;

    Point vector1 = mRightASIS - mLeftASIS;
    Point vector2 = mPubicJoin - mLeftASIS;
    Point normalRef = vector1.cross(vector2);

    mPlaneAPP.init(normalRef, mPubicJoin);

    extractFemurAxisVector(true);
	extractFemurAxisVector(false);

    isInit = true;
}

Point HipPelvis::getMidASIS() const
{
    Plane sagital;
    sagital.init(getPelvisVectorASIS(), mPubicJoin);

    Line lineASIS = Line::makeLineWithPoints(mRightASIS, mLeftASIS);
    Point midASIS = sagital.getInterceptionLinePoint(lineASIS);
    return midASIS;
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

Point HipPelvis::getFemurVectorLatMed(const Point& pCenterOfRotation) const
{
    Line myLine(mFemurAxisVector, mFemurPointInsideCenter);
    Line perpendicular = myLine.getPerpendicularLine(pCenterOfRotation);
    Point vector = perpendicular.getDirectVector();
    Plane temp;
    temp.init(vector, mFemurPointInsideCenter);
    temp.reverseByPoint(pCenterOfRotation);
    return temp.getNormalVector();
}

Point HipPelvis::getFemurVectorInfSup() const
{
    return mFemurAxisVector;
}

Point HipPelvis::getFemurVectorInfSupOppsite() const
{
	return mFemurAxisVectorOppsite;
}

double HipPelvis::getHipLengthDistance(const Point& femurHeadCenter) const
{
	Plane coronal;
	Point temp = mFemurPointInsideCenter + 10. * mFemurAxisVector;
	coronal.init(femurHeadCenter, mFemurPointInsideCenter, temp);

	Point ref;
	Line lineASIS = Line::makeLineWithPoints(coronal.getProjectionPoint(mRightASIS), coronal.getProjectionPoint(mLeftASIS));
	if (mSide == RIGHT_SIDE)
	{
		ref = coronal.getProjectionPoint(mRightLesserTrochanter);
	}
	else
	{
		ref = coronal.getProjectionPoint(mLeftLesserTrochanter);
	}

	return lineASIS.getDistanceFromPoint(ref);
}

double HipPelvis::getHipLengthDistanceOppsite(const Point& femurHeadCenterOppsite) const
{
	Plane coronal;
	Point temp = mFemurPointInsideCenterOppsite + 10. * mFemurAxisVectorOppsite;
	coronal.init(femurHeadCenterOppsite, mFemurPointInsideCenterOppsite, temp);

	Point ref;
	Line lineASIS = Line::makeLineWithPoints(coronal.getProjectionPoint(mRightASIS), coronal.getProjectionPoint(mLeftASIS));
	if (mSide == RIGHT_SIDE)
	{
		ref = coronal.getProjectionPoint(mLeftLesserTrochanter);
	}
	else
	{
		ref = coronal.getProjectionPoint(mRightLesserTrochanter);
	}

	return lineASIS.getDistanceFromPoint(ref);
}

double HipPelvis::getCombinedOffsetDistance(const Point& femurHeadCenter) const
{
	Plane coronal;
	Point temp = mFemurPointInsideCenter + 10. * mFemurAxisVector;
	coronal.init(femurHeadCenter, mFemurPointInsideCenter, temp);

	Line femoralCanalAxes(mFemurAxisVector, mFemurPointInsideCenter);
	return femoralCanalAxes.getDistanceFromPoint(coronal.getProjectionPoint(mPubicJoin));
}

double HipPelvis::getCombinedOffsetDistanceOppsite(const Point& femurHeadCenterOppsite) const
{
	Plane coronal;
	Point temp = mFemurPointInsideCenterOppsite + 10. * mFemurAxisVectorOppsite;
	coronal.init(femurHeadCenterOppsite, mFemurPointInsideCenterOppsite, temp);

	Line femoralCanalAxes(mFemurAxisVectorOppsite, mFemurPointInsideCenterOppsite);
	return femoralCanalAxes.getDistanceFromPoint(coronal.getProjectionPoint(mPubicJoin));
}

std::pair<Point, Point> HipPelvis::getAbductionAnteversionVectorsZX(const Point& pCenterOfRotation, double pAbductionAngle, double pAnteversionAngle) const
{
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
    //referenceVector = -sagital.getNormalVector(); //rotateVector will rotate towards referenceVector.
    //cv::Mat referenceVectorMat = rotMatrix * referenceVector.ToMatPoint();
    //referenceVector = Point(referenceVectorMat);
    //rotationAxis = rotateVector.cross(referenceVector); //it needs make cross because rotateVector will rotate towards referenceVector.
	
	///////////////////////////// 2nd way

	//referenceVector = getPelvisVectorInfSup(); //rotateVector will rotate towards referenceVector.
	//cv::Mat referenceVectorMat = rotMatrix * referenceVector.ToMatPoint();
	//referenceVector = Point(referenceVectorMat);
	//rotationAxis = rotateVector.cross(referenceVector);

	///////////////////////////////////////

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

double HipPelvis::getFemurVersion(const HipFemur& pFemur, const Point& pNeckAxisVectorToHead) const
{
	Point lateralMedialAxis = pFemur.getMedialEpicondyle() - pFemur.getLateralEpicondyle();
	lateralMedialAxis.normalice();
	Point midPoint = (pFemur.getMedialEpicondyle() + pFemur.getLateralEpicondyle()) / 2.;

	Plane axial;
	Point anteriorDirectionVector;

	if (pFemur.getFemurSide() == mSide)
	{
		axial.init(mFemurAxisVector, mFemurPointInsideCenter);
	}
	else
	{
		axial.init(mFemurAxisVectorOppsite, mFemurPointInsideCenterOppsite);
	}

	Point neckAxis = axial.getProjectionVector(pNeckAxisVectorToHead);
	lateralMedialAxis = axial.getProjectionVector(lateralMedialAxis);

	neckAxis.normalice();
	lateralMedialAxis.normalice();

	double angle = ImplantTools::getAngleBetweenVectors(neckAxis, lateralMedialAxis);

	if (pFemur.getFemurSide() == PelvisSide::RIGHT_SIDE)
	{
		if (pFemur.getFemurSide() == mSide)
		{
			anteriorDirectionVector = lateralMedialAxis.cross(mFemurAxisVector);
		}
		else
		{
			anteriorDirectionVector = lateralMedialAxis.cross(mFemurAxisVectorOppsite);
		}
	}
	else
	{
		if (pFemur.getFemurSide() == mSide)
		{
			anteriorDirectionVector = mFemurAxisVector.cross(lateralMedialAxis);
		}
		else
		{
			anteriorDirectionVector = mFemurAxisVectorOppsite.cross(lateralMedialAxis);
		}
	}

	anteriorDirectionVector.normalice();

	Plane coronal;
	coronal.init(anteriorDirectionVector, midPoint);

	Point checkPoint = midPoint + 1000. * neckAxis;

	if (coronal.eval(checkPoint) >= 0)
	{
		return angle;
	}
	else
	{
		return -angle;
	}
}

double HipPelvis::getFemurVersion(const HipFemur& pFemur) const
{
	Point neckAxis = getNeckAxisVector(pFemur.getHeadCenter(), pFemur.getNeck(), pFemur.getGreaterTrochanter(), pFemur.getFemur(), pFemur.getFemurSide());
	Point lateralMedialAxis = pFemur.getMedialEpicondyle() - pFemur.getLateralEpicondyle();
	lateralMedialAxis.normalice();
	Point midPoint = (pFemur.getMedialEpicondyle() + pFemur.getLateralEpicondyle()) / 2.;

	Plane axial;
	Point anteriorDirectionVector;

	if (pFemur.getFemurSide() == mSide)
	{
		axial.init(mFemurAxisVector, mFemurPointInsideCenter);
	}
	else
	{
		axial.init(mFemurAxisVectorOppsite, mFemurPointInsideCenterOppsite);
	}

	neckAxis = axial.getProjectionVector(neckAxis);
	lateralMedialAxis = axial.getProjectionVector(lateralMedialAxis);

	neckAxis.normalice();
	lateralMedialAxis.normalice();

	double angle = ImplantTools::getAngleBetweenVectors(neckAxis, lateralMedialAxis);

	if (pFemur.getFemurSide() == PelvisSide::RIGHT_SIDE)
	{
		if (pFemur.getFemurSide() == mSide)
		{
			anteriorDirectionVector = lateralMedialAxis.cross(mFemurAxisVector);
		}
		else
		{
			anteriorDirectionVector = lateralMedialAxis.cross(mFemurAxisVectorOppsite);
		}
	}
	else
	{
		if (pFemur.getFemurSide() == mSide)
		{
			anteriorDirectionVector = mFemurAxisVector.cross(lateralMedialAxis);
		}
		else
		{
			anteriorDirectionVector = mFemurAxisVectorOppsite.cross(lateralMedialAxis);
		}
	}

	anteriorDirectionVector.normalice();

	Plane coronal;
	coronal.init(anteriorDirectionVector, midPoint);
	
	Point checkPoint = midPoint + 1000. * neckAxis;

	if (coronal.eval(checkPoint) >= 0)
	{
		return angle;
	}
	else
	{
		return -angle;
	}
}

Point HipPelvis::getNeckAxisVector(const Point& femurHeadCenter, const Point& femurNeck, const Point& greaterTrochanter, const vtkSmartPointer<vtkPolyData>& femurPoly, const PelvisSide& side) const
{
	Point refLesserTrochanter, refFemurAxis, refFemurPoint;

	if (side == LEFT_SIDE)
	{
		refLesserTrochanter = mLeftLesserTrochanter;
	}
	else
	{
		refLesserTrochanter = mRightLesserTrochanter;
	}

	if (side == mSide)
	{
		refFemurAxis = mFemurAxisVector;
		refFemurPoint = mFemurPointInsideCenter;
	}
	else
	{
		refFemurAxis = mFemurAxisVectorOppsite;
		refFemurPoint = mFemurPointInsideCenterOppsite;
	}

	Point midPoint = (refLesserTrochanter + greaterTrochanter) / 2.;

	Line canalAxis(refFemurAxis, refFemurPoint);

	midPoint = canalAxis.getProjectPoint(midPoint);

	Point tempAxis = femurHeadCenter - midPoint;
	tempAxis.normalice();

	double radius = -1;
	Point neckCenter = midPoint;

	for (float i = 0.2; i < 0.9; i+= 0.1)
	{
		Point cutPoint = midPoint + i * (femurHeadCenter - midPoint);

		auto contour = ImplantTools::getContours(femurPoly, tempAxis, cutPoint);
		if (contour->GetNumberOfPoints() != 0)
		{
			std::pair<Point, double> circle = ImplantTools::minCircle(contour, tempAxis);

			if (circle.second > 0  && (circle.second < radius || radius < 0))
			{
				radius = circle.second;
				neckCenter = circle.first;
			}
		}

	}

	Point neckAxis = femurHeadCenter - neckCenter;
	neckAxis.normalice();
	return neckAxis;
}

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


vtkSmartPointer<vtkPolyData> HipPelvis::getPelvisFemurCutTest()
{
    Point vectorInfoSup = getPelvisVectorInfSup();

    Plane axial, sagital;
    Point initAxialPoint = ((mLeftLesserTrochanter + mRightLesserTrochanter) / 2.0) - 20 * vectorInfoSup;
    Point refLesserTrochanter;

    if (mSide == LEFT_SIDE)
    {
        refLesserTrochanter = mLeftLesserTrochanter;
    }
    else
    {
        refLesserTrochanter = mRightLesserTrochanter;
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
    return femurData;
}