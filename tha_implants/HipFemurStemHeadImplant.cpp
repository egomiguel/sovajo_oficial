#include "HipFemurStemHeadImplant.hpp"
#include "ImplantsException.hpp"
#include "ImplantTools.hpp"

using namespace THA::IMPLANTS;

HipFemurStemHeadImplant::HipFemurStemHeadImplant()
{
    isInit = false;
}

void HipFemurStemHeadImplant::init(const Point& pHeadBasePoint1, const Point& pHeadBasePoint2, const Point& pHeadBasePoint3, 
	const Point& pHeadInsideCenterTopPoint, const std::vector<Point>& pExternalPoints)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_FEMUR_STEM_HEAD_IMPLANT;
    }
	
	Plane headPlaneBase;
	headPlaneBase.init(pHeadBasePoint1, pHeadBasePoint2, pHeadBasePoint3);
	headPlaneBase.reverseByPoint(pHeadInsideCenterTopPoint);

	mHeadInfSupVector = headPlaneBase.getNormalVector();
	mHeadInsideCenterTopPoint = pHeadInsideCenterTopPoint;

	auto fitSphere = ImplantTools::fitSphere(pExternalPoints);
	mSphereCenter = fitSphere.first;
	mRadius = fitSphere.second;

	if (mRadius <= 0)
	{
		throw ImplantExceptionCode::CAN_NOT_FIT_SPHERE_TO_CUP_POINTS;
	}

    isInit = true;
}

Point HipFemurStemHeadImplant::getVectorInfSup() const
{
	return mHeadInfSupVector;
}

Point HipFemurStemHeadImplant::getInsideCenterTopPoint() const
{
	return mHeadInsideCenterTopPoint;
}

Point HipFemurStemHeadImplant::getCenterOfSphere() const
{
	return mSphereCenter;
}
