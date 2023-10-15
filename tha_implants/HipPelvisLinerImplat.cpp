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

void HipPelvisLinerImplant::init(const Point& pTopPoint, const std::vector<Point>& pInternalSpherePoints)
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

Point HipPelvisLinerImplant::getCenterToTopVector() const
{
	Point vector = mTopPoint - mCenter;
	vector.normalice();
	return vector;
}
