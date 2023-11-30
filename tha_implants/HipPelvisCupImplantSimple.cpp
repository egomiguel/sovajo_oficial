#include "HipPelvisCupImplantSimple.hpp"
#include "ImplantsException.hpp"
#include "ImplantTools.hpp"

using namespace THA::IMPLANTS;

HipPelvisCupImplantSimple::HipPelvisCupImplantSimple()
{
    isInit = false;
}

void HipPelvisCupImplantSimple::init(const Point& pTopPoint, const Point& pBasePoint1, const Point& pBasePoint2, const Point& pBasePoint3, const std::vector<Point>& pExternalPoints)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_PELVIS_CUT_IMPLANT;
    }

	Plane mBasePlane;
	mBasePlane.init(pBasePoint1, pBasePoint2, pBasePoint3);
	mBasePlane.reverseByPoint(pTopPoint);
	mVectorZ = mBasePlane.getNormalVector();

	if (pExternalPoints.size() > 0)
	{
		auto fitSphere = ImplantTools::fitSphere(pExternalPoints);
		mCenter = fitSphere.first;
	}
	else
	{
		mCenter = mBasePlane.getProjectionPoint(pTopPoint);
	}

    isInit = true;
}

void HipPelvisCupImplantSimple::init(const Point& pVectorZ, const std::vector<Point>& pExternalPoints)
{
	if (isInit == true)
	{
		throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_PELVIS_CUT_IMPLANT;
	}

	mVectorZ = pVectorZ;

	if (pExternalPoints.size() > 0)
	{
		auto fitSphere = ImplantTools::fitSphere(pExternalPoints);
		mCenter = fitSphere.first;
	}
	else
	{
		mCenter = Point(0, 0, 0);
	}

	isInit = true;
}

void HipPelvisCupImplantSimple::init(const Point& pVectorZ, const Point& pCenter)
{
	if (isInit == true)
	{
		throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_PELVIS_CUT_IMPLANT;
	}

	mVectorZ = pVectorZ;
	mCenter = pCenter;

	isInit = true;
}

Point HipPelvisCupImplantSimple::getVectorZ() const
{
    return mVectorZ;
}

Point HipPelvisCupImplantSimple::getCenterOfRotationImplant() const
{
	return mCenter;
}