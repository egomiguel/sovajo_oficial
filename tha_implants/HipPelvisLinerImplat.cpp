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

void HipPelvisLinerImplant::init(const Point& pTopPoint, const Point& pBasePoint1, const Point& pBasePoint2, const Point& pBasePoint3, const std::vector<Point>& pInternalSpherePoints)
{
	if (isInit == true)
	{
		throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_PELVIS_LINER_IMPLANT;
	}

	this->mTopPoint = pTopPoint;
	this->mBasePoint1 = pBasePoint1;
	this->mBasePoint2 = pBasePoint2;
	this->mBasePoint3 = pBasePoint3;

	auto fitSphere = ImplantTools::fitSphere(pInternalSpherePoints);
	mCenter = fitSphere.first;
	mRadius = fitSphere.second;

	if (mRadius <= 0)
	{
		throw ImplantExceptionCode::CAN_NOT_FIT_SPHERE_TO_CUP_POINTS;
	}

	mBasePlane.init(pBasePoint1, pBasePoint2, pBasePoint3);
	mBasePlane.reverseByPoint(pTopPoint);

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

Plane HipPelvisLinerImplant::getBasePlane() const
{
	return mBasePlane;
}
