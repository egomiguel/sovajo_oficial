
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
#include <iostream>
#include "FemurImplantOnePlane.hpp"
#include "ImplantsException.hpp"


using namespace UKA::IMPLANTS;

FemurImplantOnePlane::FemurImplantOnePlane()
{
}

void FemurImplantOnePlane::init(const Plane& pPosterior, const Point& pRodBasePoint, const Point& pRodTopPoint,
						const std::vector<Point>& pSortPointsSide1, const std::vector<Point>& pSortPointsSide2,
						const vtkSmartPointer<vtkPolyData> implantModel, const FemurImplantInfo& pImplantInfo)
{
	FemurImplant::init(pPosterior, implantModel, pImplantInfo);

	if (pSortPointsSide1.size() < 3 || pSortPointsSide2.size() < 3)
	{
		throw ImplantExceptionCode::FEMUR_IMPLANT_MUST_HAVE_ENOUGH_SIDE_POINTS;
	}

    this->mRodBasePoint = pRodBasePoint;
    this->mRodTopPoint = pRodTopPoint;
    this->mSideBorder1 = pSortPointsSide1;
	this->mSideBorder2 = pSortPointsSide2;

	mPosterior.reverseByPoint(mRodTopPoint);

	if (mPosterior.getDistanceFromPoint(*mSideBorder1.begin()) > mPosterior.getDistanceFromPoint(*mSideBorder1.rbegin()))
	{
		std::reverse(mSideBorder1.begin(), mSideBorder1.end());
	}

	if (mPosterior.getDistanceFromPoint(*mSideBorder2.begin()) > mPosterior.getDistanceFromPoint(*mSideBorder2.rbegin()))
	{
		std::reverse(mSideBorder2.begin(), mSideBorder2.end());
	}

	bool tResult;
	Plane side1 = Plane::getBestPlane(pSortPointsSide1, tResult);
	Plane side2 = Plane::getBestPlane(pSortPointsSide2, tResult);
	/*
		Both vectors must have the same direction to calculate an average vector with best precision.
	*/
	side1.reverseByPoint(mRodBasePoint);
	side2.reverseByPoint(mRodBasePoint, false);
	mSizeMidVector = side1.getNormalVector() + side2.getNormalVector();
	mSizeMidVector.normalice();

	Point crossBase = mPosterior.getNormalVector().cross(mSizeMidVector);
	mDistal.init(crossBase, mRodBasePoint);
	mDistal.reverseByPoint(mRodTopPoint);

	mVectorAP = mPosterior.getNormalVector();
	mVectorForceLine = mDistal.getNormalVector();
	mVectorTEA = mVectorForceLine.cross(mVectorAP);
	mVectorTEA.normalice();

	mRodTopPointProjectedOnBase = mDistal.getProjectionPoint(pRodTopPoint);
}

FemurImplantOnePlane::FemurImplantOnePlane(const FemurImplantOnePlane& pImplant)
{
	this->mPosterior = pImplant.mPosterior;
	this->mRodBasePoint = pImplant.mRodBasePoint;
	this->mRodTopPoint = pImplant.mRodTopPoint;
	this->mSideBorder1 = pImplant.mSideBorder1;
	this->mSideBorder2 = pImplant.mSideBorder2;
    this->mImplantModel = pImplant.mImplantModel;
    this->isInit = pImplant.isInit;
    this->mImplantInfo = pImplant.mImplantInfo;
	this->mVectorAP = pImplant.mVectorAP;
	this->mVectorForceLine = pImplant.mVectorForceLine;
	this->mVectorTEA = pImplant.mVectorTEA;
	this->mSizeMidVector = pImplant.mSizeMidVector;
	this->mDistal = pImplant.mDistal;
	this->mRodTopPointProjectedOnBase = pImplant.mRodTopPointProjectedOnBase;
}

Plane FemurImplantOnePlane::getDistalPlane() const
{
	return mDistal;
}

Plane FemurImplantOnePlane::getMidPlane() const
{
	auto const count1 = static_cast<float>(mSideBorder1.size());
	auto const count2 = static_cast<float>(mSideBorder2.size());

	Point sum1, sum2;

	auto it1 = mSideBorder1.begin();
	auto it2 = mSideBorder1.end();

	for ( ; it1 != it2; ++it1 )
	{
		sum1 += *it1;
	}

	it1 = mSideBorder2.begin();
	it2 = mSideBorder2.end();

	for (; it1 != it2; ++it1)
	{
		sum2 += *it1;
	}

	sum1 = sum1 / count1;
	sum2 = sum2 / count2;

	Plane midPlane;
	midPlane.init(mSizeMidVector, sum1);
	midPlane.reverseByPoint(sum2);
	double distance = midPlane.getDistanceFromPoint(sum2) / 2.;
	midPlane.movePlaneOnNormal(distance);
	return midPlane;
}

double FemurImplantOnePlane::getWidthSize() const
{
	auto const count1 = static_cast<float>(mSideBorder1.size());
	auto const count2 = static_cast<float>(mSideBorder2.size());

	Point sum1, sum2;

	auto it1 = mSideBorder1.begin();
	auto it2 = mSideBorder1.end();

	for (; it1 != it2; ++it1)
	{
		sum1 += *it1;
	}

	it1 = mSideBorder2.begin();
	it2 = mSideBorder2.end();

	for (; it1 != it2; ++it1)
	{
		sum2 += *it1;
	}

	sum1 = sum1 / count1;
	sum2 = sum2 / count2;

	Plane temp;
	temp.init(mSizeMidVector, sum1);
	double distance = temp.getDistanceFromPoint(sum2);
	return distance;
}

Point FemurImplantOnePlane::getRodTopPoint() const
{
	return mRodTopPoint;
}

Point FemurImplantOnePlane::getRodBasePoint() const
{
	return mRodBasePoint;
}

Point FemurImplantOnePlane::getRodTopPointProjectedOnBase() const
{
	return mRodTopPointProjectedOnBase;
}

Point FemurImplantOnePlane::getRodTopPointProjectedOnBaseExterior() const
{
	return mRodTopPointProjectedOnBase - mImplantInfo.femurDistalThickness * mDistal.getNormalVector();
}

std::vector<Point> FemurImplantOnePlane::getSortPointsSide1() const
{
	return mSideBorder1;
}

std::vector<Point> FemurImplantOnePlane::getSortPointsSide2() const
{
	return mSideBorder2;
}

std::vector<Point> FemurImplantOnePlane::getAllSidePointsInOrder(double pOffset) const
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

std::vector<Point> FemurImplantOnePlane::getAllSidePointsInOrder(cv::Mat& pRotation, cv::Mat& pTranslation, double pOffset) const
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

Plane FemurImplantOnePlane::getBestPlaneToCurvePoints() const
{
	std::vector<Point> points = getAllSidePointsInOrder();
	bool error;
	return Plane::getBestPlane(points, error);
}

Point FemurImplantOnePlane::getDirectVectorFemurAxis() const
{
	return mVectorForceLine;
}

Point FemurImplantOnePlane::getDirectVectorTEA() const
{
	return mVectorTEA;
}

Point FemurImplantOnePlane::getDirectVectorAP() const
{
	return mVectorAP;
}
