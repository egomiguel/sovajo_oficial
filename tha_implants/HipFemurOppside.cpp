#include "HipFemurOppside.hpp"
#include "Plane.hpp"
#include "ImplantTools.hpp"
#include "ImplantsException.hpp"

using namespace THA::IMPLANTS;

HipFemurOppside::HipFemurOppside()
{
	isInit = false;
}

void HipFemurOppside::init(const Point& headCenter, const Point& greaterTrochanter, const Point& lesserTrochanter,
							const Point& femurKneeCenter)
{
	if (isInit == true)
	{
		throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_FEMUR;
	}

	mHeadCenter = headCenter;
	mGreaterTrochanter = greaterTrochanter;
	mLesserTrochanter = lesserTrochanter;
	mCanalAxisVectorInfSup = greaterTrochanter - femurKneeCenter;
	mCanalAxisVectorInfSup.normalice();
	mCanalAxisPoint = femurKneeCenter;
	mKneeCenter = femurKneeCenter;
	isInit = true;
}

Point HipFemurOppside::getHeadCenter() const
{
	return mHeadCenter;
}

Point HipFemurOppside::getGreaterTrochanter() const
{
	return mGreaterTrochanter;
}

Point HipFemurOppside::getLesserTrochanter() const
{
	return mLesserTrochanter;
}

Point HipFemurOppside::getCanalAxisVectorInfSup() const
{
	return mCanalAxisVectorInfSup;
}

Point HipFemurOppside::getCanalAxisPoint() const
{
	return mCanalAxisPoint;
}

Point HipFemurOppside::getKneeCenter() const
{
	return mKneeCenter;
}
