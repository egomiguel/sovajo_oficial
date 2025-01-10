#include "HipPelvisLinerImplat.hpp"
#include "ImplantsException.hpp"
#include "ImplantTools.hpp"
#include "vtkSphereSource.h"
#include "vtkPlane.h"
#include "vtkPlaneCollection.h"
#include "vtkClipClosedSurface.h"

using namespace THA::IMPLANTS;

HipPelvisLinerImplant::HipPelvisLinerImplant()
{
	isInit = false;
}

void HipPelvisLinerImplant::init(const Point& pTopPoint, const std::vector<Point>& pInternalSpherePoints, double pThickness)
{
	if (isInit == true)
	{
		throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_PELVIS_LINER_IMPLANT;
	}

	this->mTopPoint = pTopPoint;

	auto fitSphere = ImplantTools::fitSphere(pInternalSpherePoints);
	mCenter = fitSphere.first;
	mRadius = fitSphere.second;

	if (mRadius <= 0)
	{
		throw ImplantExceptionCode::CAN_NOT_FIT_SPHERE_TO_CUP_POINTS;
	}

	mRadius += pThickness;

	isInit = true;
}

void HipPelvisLinerImplant::init(const Point& pTopPoint, const Point& pInternalCenter, double pRadius, double pThickness)
{
	if (isInit == true)
	{
		throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_PELVIS_LINER_IMPLANT;
	}

	this->mTopPoint = pTopPoint;

	mCenter = pInternalCenter;
	mRadius = pRadius;

	mRadius += pThickness;

	isInit = true;
}

Point HipPelvisLinerImplant::getTopPoint() const
{
	return mTopPoint;
}

Point HipPelvisLinerImplant::getCenterOfRotationImplant() const
{
	return mCenter;
}

double HipPelvisLinerImplant::getExternalRadius() const
{
	return mRadius;
}

Point HipPelvisLinerImplant::getCenterToTopVector() const
{
	Point vector = mTopPoint - mCenter;
	vector.normalice();
	return vector;
}
