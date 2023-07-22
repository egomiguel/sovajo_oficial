#include "HipFemurStemHeadImplant.hpp"
#include "ImplantsException.hpp"

using namespace THA::IMPLANTS;

HipFemurStemHeadImplant::HipFemurStemHeadImplant()
{
    isInit = false;
}

void HipFemurStemHeadImplant::init(const Point& pHeadBasePoint1, const Point& pHeadBasePoint2, const Point& pHeadBasePoint3, const Point& pHeadInsideCenterTopPoint)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_FEMUR_STEM_HEAD_IMPLANT;
    }
	
	Plane headPlaneBase;
	headPlaneBase.init(pHeadBasePoint1, pHeadBasePoint2, pHeadBasePoint3);

	mHeadInfSupVector = pHeadInsideCenterTopPoint - headPlaneBase.getProjectionPoint(pHeadInsideCenterTopPoint);
	mHeadInfSupVector.normalice();
	mHeadInsideCenterTopPoint = pHeadInsideCenterTopPoint;

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
