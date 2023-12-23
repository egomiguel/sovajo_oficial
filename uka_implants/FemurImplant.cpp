
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
#include <iostream>
#include "FemurImplant.hpp"
#include "ImplantsException.hpp"


using namespace UKA::IMPLANTS;

FemurImplant::FemurImplant()
{
	isInit = false;
}

void FemurImplant::init(const Plane& pPosterior, const Point& pRodBasePoint, const Point& pRodTopPoint,
						const std::vector<Point>& pSortPointsSide1, const std::vector<Point>& pSortPointsSide2,
						const vtkSmartPointer<vtkPolyData> implantModel, const FemurImplantInfo& pImplantInfo)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_FEMUR_IMPLANT;
    }

	if (pSortPointsSide1.size() < 3 || pSortPointsSide2.size() < 3)
	{
		throw ImplantExceptionCode::FEMUR_IMPLANT_MUST_HAVE_ENOUGH_SIDE_POINTS;
	}

    this->mPosterior = pPosterior;
    this->mRodBasePoint = pRodBasePoint;
    this->mRodTopPoint = pRodTopPoint;
    this->mSideBorder1 = pSortPointsSide1;
	this->mSideBorder2 = pSortPointsSide2;
    this->mImplantModel = implantModel;
    this->mImplantInfo = pImplantInfo;

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

    isInit = true;
}

FemurImplant::FemurImplant(const FemurImplant& pImplant)
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
}

FemurImplantInfo FemurImplant::getImplantInfo() const
{
    return mImplantInfo;
}

Plane FemurImplant::getPosterior() const
{
	return mPosterior;
}

Plane FemurImplant::getDistalPlane() const
{
	return mDistal;
}

Plane FemurImplant::getMidPlane() const
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

double FemurImplant::getWidthSize() const
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

Point FemurImplant::getRodTopPoint() const
{
	return mRodTopPoint;
}

Point FemurImplant::getRodBasePoint() const
{
	return mRodBasePoint;
}

std::vector<Point> FemurImplant::getSortPointsSide1() const
{
	return mSideBorder1;
}

std::vector<Point> FemurImplant::getSortPointsSide2() const
{
	return mSideBorder2;
}

std::vector<Point> FemurImplant::getAllSidePointsInOrder(double pOffset) const
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

std::vector<Point> FemurImplant::getAllSidePointsInOrder(cv::Mat& pRotation, cv::Mat& pTranslation, double pOffset) const
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

Plane FemurImplant::getBestPlaneToCurvePoints() const
{
	std::vector<Point> points = getAllSidePointsInOrder();
	bool error;
	return Plane::getBestPlane(points, error);
}

Point FemurImplant::getDirectVectorFemurAxis() const
{
	return mVectorForceLine;
}

Point FemurImplant::getDirectVectorTEA() const
{
	return mVectorTEA;
}

Point FemurImplant::getDirectVectorAP() const
{
	return mVectorAP;
}

vtkSmartPointer<vtkPolyData> FemurImplant::GetImplantModel() const
{
    return mImplantModel;
}

vtkSmartPointer<vtkPolyData> FemurImplant::GetTransformImplantModel(const itk::Rigid3DTransform<>::Pointer transform) const
{
    itk::Matrix< double, 3, 3 > rotation = transform->GetMatrix();
    itk::Vector< double, 3 > translate = transform->GetOffset();

    vtkNew<vtkMatrix4x4> m;
    m->SetElement(0, 0, rotation[0][0]);
    m->SetElement(0, 1, rotation[0][1]);
    m->SetElement(0, 2, rotation[0][2]);

    m->SetElement(1, 0, rotation[1][0]);
    m->SetElement(1, 1, rotation[1][1]);
    m->SetElement(1, 2, rotation[1][2]);

    m->SetElement(2, 0, rotation[2][0]);
    m->SetElement(2, 1, rotation[2][1]);
    m->SetElement(2, 2, rotation[2][2]);

    m->SetElement(3, 0, 0);
    m->SetElement(3, 1, 0);
    m->SetElement(3, 2, 0);
    m->SetElement(3, 3, 1);

    m->SetElement(0, 3, translate[0]);
    m->SetElement(1, 3, translate[1]);
    m->SetElement(2, 3, translate[2]);

    vtkNew<vtkTransform> vtkTransform;
    vtkTransform->SetMatrix(m);

    vtkNew<vtkTransformPolyDataFilter> transformFilter;
    transformFilter->SetInputData(mImplantModel);
    transformFilter->SetTransform(vtkTransform);
    transformFilter->Update();
    return transformFilter->GetOutput();
}