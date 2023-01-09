#include "HipFemurStemImplant.hpp"
#include "ImplantsException.hpp"


HipFemurStemImplant::HipFemurStemImplant()
{
    isInit = false;
}

void HipFemurStemImplant::init(const Point& pTopPoint, const Point& pBasePoint, const Point& pHeadCenter)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_FEMUR_STEM_IMPLANT;
    }

    this->mTopPoint = pTopPoint;
    this->mBasePoint = pBasePoint;
    this->mHeadCenter = pHeadCenter;
    isInit = true;
}

Point HipFemurStemImplant::getVectorInfoSup() const
{
    Point vector = mTopPoint - mBasePoint;
    vector.normalice();
    return vector;
}

Point HipFemurStemImplant::getVectorLatMed() const
{
    Line myLine(getVectorInfoSup(), mTopPoint);
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