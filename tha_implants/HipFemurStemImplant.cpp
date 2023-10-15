#include "HipFemurStemImplant.hpp"
#include "ImplantsException.hpp"

using namespace THA::IMPLANTS;

HipFemurStemImplant::HipFemurStemImplant()
{
    isInit = false;
}

void HipFemurStemImplant::init(const Point& pTopPoint, const Point& pBasePoint, const Point& pHeadCenter, const std::vector<Point>& pHeadPoints)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_FEMUR_STEM_IMPLANT;
    }

    this->mTopPoint = pTopPoint;
    this->mBasePoint = pBasePoint;
    this->mHeadCenter = pHeadCenter;
	bool result;
	this->mStemHeadPlane = Plane::getBestPlane(pHeadPoints, result);
	if (result == false)
	{
		throw ImplantExceptionCode::CAN_NOT_FIT_PLANE_TO_STEM_HEAD;
	}
	this->mStemHeadPlane.reverseByPoint(mTopPoint, false);
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
	return mStemHeadPlane.getNormalVector();
}

Point HipFemurStemImplant::getVectorNeckToHeadPerpendicularToInfSup() const
{
	Plane base;
	base.init(getVectorInfSup(), mTopPoint);
	Point result = base.getProjectionVector(getVectorNeckToHead());
	result.normalice();
	return result;
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

Point HipFemurStemImplant::getCanalAxisTopPoint() const
{
	return mTopPoint;
}