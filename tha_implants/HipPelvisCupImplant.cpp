#include "HipPelvisCupImplant.hpp"
#include "ImplantsException.hpp"


using namespace THA::IMPLANTS;

HipPelvisCupImplant::HipPelvisCupImplant()
{
    isInit = false;
}

void HipPelvisCupImplant::init(const Point& pTopPoint, const Point& pBasePoint1, const Point& pBasePoint2, const Point& pBasePoint3)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_PELVIS_CUT_IMPLANT;
    }

    this->mTopPoint = pTopPoint;
    this->mBasePoint1 = pBasePoint1;
    this->mBasePoint2 = pBasePoint2;
    this->mBasePoint3 = pBasePoint3;

    isInit = true;
}

Point HipPelvisCupImplant::getVectorX() const
{
    Point vector = mBasePoint1 - mBasePoint2;
    vector.normalice();
    return vector;
}

Point HipPelvisCupImplant::getVectorZ() const
{
    Plane myPlane;
    myPlane.init(mBasePoint1, mBasePoint2, mBasePoint3);
    myPlane.reverseByPoint(mTopPoint);
    return myPlane.getNormalVector();
}

Point HipPelvisCupImplant::getTopPoint() const
{
    return mTopPoint;
}

Point HipPelvisCupImplant::getCenterOfRotationImplant() const
{
    Plane myPlane;
    myPlane.init(mBasePoint1, mBasePoint2, mBasePoint3);
    return myPlane.getProjectionPoint(mTopPoint);
}

Plane HipPelvisCupImplant::getBasePlane() const
{
    Plane myPlane;
    myPlane.init(mBasePoint1, mBasePoint2, mBasePoint3);
    return myPlane;
}