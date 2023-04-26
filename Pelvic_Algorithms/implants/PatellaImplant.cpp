
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
#include <iostream>
#include "PatellaImplant.hpp"
#include "ImplantsException.hpp"

PatellaImplant::PatellaImplant()
{
	isInit = false;
}

void PatellaImplant::init(const Point& basePoint1, const Point& basePoint2, const Point& basePoint3, const Point& topCentralPoint)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_PATELLA_IMPLANT;
    }

    this->mBasePoint1 = basePoint1;
    this->mBasePoint2 = basePoint2;
    this->mBasePoint3 = basePoint3;
    this->mTopCentralPoint = topCentralPoint;

    mBasePlane.init(mBasePoint1, mBasePoint2, mBasePoint3);
    mBasePlane.reverseByPoint(mTopCentralPoint);

    isInit = true;
}

PatellaImplant::PatellaImplant(const PatellaImplant& pImplant)
{
    this->mBasePlane = pImplant.mBasePlane;
    this->mBasePoint1 = pImplant.mBasePoint1;
    this->mBasePoint2 = pImplant.mBasePoint2;
    this->isInit = pImplant.isInit;
    this->mBasePoint3 = pImplant.mBasePoint3;
    this->mTopCentralPoint = pImplant.mTopCentralPoint;
}

Plane PatellaImplant::getBasePlane() const
{
    return mBasePlane;
}

Point PatellaImplant::getNormalVector() const
{
    return mBasePlane.getNormalVector();
}

Point PatellaImplant::getCentralPointOnBase() const
{
    return mBasePlane.getProjectionPoint(mTopCentralPoint);
}

double PatellaImplant::getThickness() const
{
    return mBasePlane.getDistanceFromPoint(mTopCentralPoint);
}

Point PatellaImplant::getTopPoint() const
{
    return mTopCentralPoint;
}



