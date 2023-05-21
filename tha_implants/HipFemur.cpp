#include "HipFemur.hpp"
#include "Plane.hpp"
#include "ImplantTools.hpp"
#include "ImplantsException.hpp"

using namespace THA::IMPLANTS;

HipFemur::HipFemur()
{
	isInit = false;
}

void HipFemur::init(const Point& headCenter, const Point& neck, const Point& greaterTrochanter, const Point& lesserTrochanter,
	const Point& medialEpicondyle, const Point& lateralEpicondyle, const vtkSmartPointer<vtkPolyData>& femurPoly,
	const PelvisSide& femurSide)
{
	if (isInit == true)
	{
		throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_FEMUR;
	}

	mHeadCenter = headCenter;
	mNeck = neck;
	mGreaterTrochanter = greaterTrochanter;
	mLesserTrochanter = lesserTrochanter;
	mMedialEpicondyle = medialEpicondyle;
	mLateralEpicondyle = lateralEpicondyle;
	mFemur = femurPoly;
	mFemurSide = femurSide;
	isInit = true;
}

vtkSmartPointer<vtkPolyData> HipFemur::getFemur() const
{
	return mFemur;
}

Point HipFemur::getHeadCenter() const
{
	return mHeadCenter;
}

Point HipFemur::getNeck() const
{
	return mNeck;
}

Point HipFemur::getGreaterTrochanter() const
{
	return mGreaterTrochanter;
}

Point HipFemur::getLesserTrochanter() const
{
	return mLesserTrochanter;
}

Point HipFemur::getMedialEpicondyle() const
{
	return mMedialEpicondyle;
}

Point HipFemur::getLateralEpicondyle() const
{
	return mLateralEpicondyle;
}

PelvisSide HipFemur::getFemurSide() const
{
	return mFemurSide;
}
