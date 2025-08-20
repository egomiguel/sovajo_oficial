
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
#include <iostream>
#include "FemurImplantThreePlane.hpp"
#include "ImplantsException.hpp"
#include "ImplantTools.hpp"

using namespace UKA::IMPLANTS;

FemurImplantThreePlane::FemurImplantThreePlane()
{
}

void FemurImplantThreePlane::init(const Plane& pPosterior, const Plane& pCenter, const Plane& pAnterior, const Point& pRodTopPoint,
	const Point& pRodBaseExtremeSide1, const Point& pRodBaseExtremeSide2,
	const vtkSmartPointer<vtkPolyData> implantModel, const FemurImplantInfo& pImplantInfo)
{
	FemurImplant::init(pPosterior, implantModel, pImplantInfo);

	this->mCenter = pCenter;
	this->mAnterior = pAnterior;
    this->mRodTopPoint = pRodTopPoint;
    this->mRodBasePoint1 = pRodBaseExtremeSide1;
	this->mRodBasePoint2 = pRodBaseExtremeSide2;

	Point midPoint = (mRodBasePoint1 + mRodBasePoint2) / 2;
	mPosterior.reverseByPoint(mRodTopPoint);
	mCenter.reverseByPoint(mRodTopPoint);
	mAnterior.reverseByPoint(mRodTopPoint);
	mVectorAP = mPosterior.getNormalVector();
	mVectorTEA = mCenter.getNormalVector().cross(mPosterior.getNormalVector());
	mVectorTEA.normalice();
	Point crossBase = mVectorAP.cross(mVectorTEA);
	mDistal.init(crossBase, midPoint);
	mDistal.reverseByPoint(mRodTopPoint);
	mVectorForceLine = mDistal.getNormalVector();
}

FemurImplantThreePlane::FemurImplantThreePlane(const FemurImplantThreePlane& pImplant)
{
	this->mPosterior = pImplant.mPosterior;
	this->mRodTopPoint = pImplant.mRodTopPoint;
	this->mRodBasePoint1 = pImplant.mRodBasePoint1;
	this->mRodBasePoint2 = pImplant.mRodBasePoint2;
    this->mImplantModel = pImplant.mImplantModel;
    this->isInit = pImplant.isInit;
    this->mImplantInfo = pImplant.mImplantInfo;
	this->mVectorAP = pImplant.mVectorAP;
	this->mVectorForceLine = pImplant.mVectorForceLine;
	this->mVectorTEA = pImplant.mVectorTEA;
	this->mDistal = pImplant.mDistal;
	this->mCenter = pImplant.mCenter;
	this->mAnterior = pImplant.mAnterior;
}

Plane FemurImplantThreePlane::getDistalPlane() const
{
	return mDistal;
}

Plane FemurImplantThreePlane::getMidPlane() const
{
	Point midPoint = (mRodBasePoint1 + mRodBasePoint2) / 2;
	Plane midPlane;
	midPlane.init(mVectorTEA, midPoint);
	return midPlane;
}

Point FemurImplantThreePlane::getRotationPoint() const
{
	return (mRodBasePoint1 + mRodBasePoint2) / 2;
}

double FemurImplantThreePlane::getWidthSize() const
{
	return ImplantTools::getDistanceBetweenPoints(mRodBasePoint1, mRodBasePoint2);
}

Point FemurImplantThreePlane::getRodTopPoint() const
{
	return mRodTopPoint;
}

Point FemurImplantThreePlane::getRodTopPointProjectedOnBase() const
{
	Line lineTemp = Line::makeLineWithPoints(mRodBasePoint1, mRodBasePoint2);
	return lineTemp.getProjectPoint(mRodTopPoint);
}

Point FemurImplantThreePlane::getRodTopPointProjectedOnBaseExterior() const
{
	return getRodTopPointProjectedOnBase() - mImplantInfo.femurDistalThickness * mDistal.getNormalVector();
}

/*
std::vector<Point> FemurImplantThreePlane::getSortPointsSide1() const
{
	return mSideBorder1;
}

std::vector<Point> FemurImplantThreePlane::getSortPointsSide2() const
{
	return mSideBorder2;
}

std::vector<Point> FemurImplantThreePlane::getAllSidePointsInOrder(double pOffset) const
{
	std::vector<Point> result;

	Point move1, move2;

	if (pOffset != 0)
	{
		bool error;
		Plane sidePlane1 = Plane::getBestPlane(mSideBorder1, error);
		sidePlane1.reverseByPoint(mRodBasePoint, false);
		move1 = sidePlane1.getNormalVector();

		Plane sidePlane2 = Plane::getBestPlane(mSideBorder2, error);
		sidePlane2.reverseByPoint(mRodBasePoint, false);
		move2 = sidePlane2.getNormalVector();
	}

	for (auto it = mSideBorder1.rbegin(); it != mSideBorder1.rend(); ++it)
	{
		Point temp = (*it) + pOffset * move1;
		result.push_back(temp);
	}

	for (auto& it : mSideBorder2)
	{
		Point temp = it + pOffset * move2;
		result.push_back(temp);
	}

	return result;
}

std::vector<Point> FemurImplantThreePlane::getAllSidePointsInOrder(cv::Mat& pRotation, cv::Mat& pTranslation, double pOffset) const
{
	std::vector<Point> result;

	Point move1, move2;

	if (pOffset != 0)
	{
		bool error;
		Plane sidePlane1 = Plane::getBestPlane(mSideBorder1, error);
		sidePlane1.reverseByPoint(mRodBasePoint, false);
		move1 = sidePlane1.getNormalVector();

		Plane sidePlane2 = Plane::getBestPlane(mSideBorder2, error);
		sidePlane2.reverseByPoint(mRodBasePoint, false);
		move2 = sidePlane2.getNormalVector();
	}

	for (auto it = mSideBorder1.rbegin(); it != mSideBorder1.rend(); ++it)
	{
		Point temp = (*it) + pOffset * move1;
		cv::Mat tempMat = pRotation * temp.ToMatPoint() + pTranslation;
		result.push_back(Point(tempMat));
	}

	for (auto& it : mSideBorder2)
	{
		Point temp = it + pOffset * move2;
		cv::Mat tempMat = pRotation * temp.ToMatPoint() + pTranslation;
		result.push_back(Point(tempMat));
	}

	return result;
}

Plane FemurImplantThreePlane::getBestPlaneToCurvePoints() const
{
	std::vector<Point> points = getAllSidePointsInOrder();
	bool error;
	return Plane::getBestPlane(points, error);
}

*/

Point FemurImplantThreePlane::getDirectVectorFemurAxis() const
{
	return mVectorForceLine;
}

Point FemurImplantThreePlane::getDirectVectorTEA() const
{
	return mVectorTEA;
}

Point FemurImplantThreePlane::getDirectVectorAP() const
{
	return mVectorAP;
}

Plane FemurImplantThreePlane::getAnteriorPlane() const
{
	return mAnterior;
}

Plane FemurImplantThreePlane::getCenterPlane() const
{
	return mCenter;
}
