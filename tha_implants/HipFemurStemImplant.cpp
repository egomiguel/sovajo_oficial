#include "HipFemurStemImplant.hpp"
#include "ImplantsException.hpp"

using namespace THA::IMPLANTS;

HipFemurStemImplant::HipFemurStemImplant()
{
    isInit = false;
}

void HipFemurStemImplant::init(const Point& pRodCenter, const Point& pBasePoint, const Point& pHeadCenter, const std::vector<Point>& pHeadPoints)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_FEMUR_STEM_IMPLANT;
    }

    this->mRodCenter = pRodCenter;
    this->mBasePoint = pBasePoint;
    this->mHeadCenter = pHeadCenter;
	bool result;
	this->mStemHeadPlane = Plane::getBestPlane(pHeadPoints, result);
	if (result == false)
	{
		throw ImplantExceptionCode::CAN_NOT_FIT_PLANE_TO_STEM_HEAD;
	}
	this->mStemHeadPlane.reverseByPoint(pRodCenter, false);
    isInit = true;
}

Point HipFemurStemImplant::getBasePoint() const
{
	return mBasePoint;
}

Point HipFemurStemImplant::getVectorInfSup() const
{
    Point vector = mRodCenter - mBasePoint;
    vector.normalice();
    return vector;
}

Point HipFemurStemImplant::getVectorNeckToHead() const
{
	return mStemHeadPlane.getNormalVector();
}

Point HipFemurStemImplant::getVectorNeckToHeadPerpendicularToInfSup() const
{
	Plane base;
	base.init(getVectorInfSup(), mRodCenter);
	Point result = base.getProjectionVector(getVectorNeckToHead());
	result.normalice();
	return result;
}

Point HipFemurStemImplant::getVectorLatMed() const
{
    Line myLine(getVectorInfSup(), mRodCenter);
    Line perpendicular = myLine.getPerpendicularLine(mHeadCenter);
    Point vector = perpendicular.getDirectVector();
    Plane temp;
    temp.init(vector, mRodCenter);
    temp.reverseByPoint(mHeadCenter);
    return temp.getNormalVector();
}

Point HipFemurStemImplant::getHeadCenter() const
{
    return mHeadCenter;
}

Point HipFemurStemImplant::getCanalAxisRodCenter() const
{
	return mRodCenter;
}