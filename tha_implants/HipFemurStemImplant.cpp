#include "HipFemurStemImplant.hpp"
#include "ImplantsException.hpp"

using namespace THA::IMPLANTS;

HipFemurStemImplant::HipFemurStemImplant()
{
    isInit = false;
}

void HipFemurStemImplant::init(const Point& pTopPoint, const Point& pBasePoint, const Point& pHeadCenter, const Point& pNeckCenter)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_FEMUR_STEM_IMPLANT;
    }

    this->mTopPoint = pTopPoint;
    this->mBasePoint = pBasePoint;
    this->mHeadCenter = pHeadCenter;
	this->mNeckCenter = pNeckCenter;
    isInit = true;
}

Point HipFemurStemImplant::getVectorInfSup() const
{
    Point vector = mTopPoint - mBasePoint;
    vector.normalice();
    return vector;
}

Point HipFemurStemImplant::getVectorNeckToHead() const
{
	Point neckAxis = mHeadCenter - mNeckCenter;
	neckAxis.normalice();
	return neckAxis;
}

Point HipFemurStemImplant::getVectorLatMed() const
{
    Line myLine(getVectorInfSup(), mTopPoint);
    Line perpendicular = myLine.getPerpendicularLine(mHeadCenter);
    Point vector = perpendicular.getDirectVector();
    Plane temp;
    temp.init(vector, mTopPoint);
    temp.reverseByPoint(mHeadCenter);
    return temp.getNormalVector();
}

Point HipFemurStemImplant::getHeadCenter() const
{
    return mHeadCenter;
}

Point HipFemurStemImplant::getCanalAxisPoint() const
{
	return mBasePoint;
}