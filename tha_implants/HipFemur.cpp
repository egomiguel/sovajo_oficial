#include "HipFemur.hpp"
#include "ImplantTools.hpp"
#include "ImplantsException.hpp"

using namespace THA::IMPLANTS;

HipFemur::HipFemur()
{
	isInit = false;
}

void HipFemur::init(const Point& headCenter, const Point& neckCenter, const Point& canalCenter, const Point& lesserTrochanter,
	const Point& medialEpicondyle, const Point& lateralEpicondyle, const Point& femurKneeCenter, 
	const vtkSmartPointer<vtkPolyData>& femurPoly)
{
	if (isInit == true)
	{
		throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_FEMUR;
	}

	mHeadCenter = headCenter;
	mNeckCenter = neckCenter;
	//mGreaterTrochanter = greaterTrochanter;
	mLesserTrochanter = lesserTrochanter;
	mMedialEpicondyle = medialEpicondyle;
	mLateralEpicondyle = lateralEpicondyle;
	mFemur = femurPoly;
	mCanalAxisPoint = canalCenter;
	mCanalAxisVectorInfSup = mCanalAxisPoint - femurKneeCenter;
	mCanalAxisVectorInfSup.normalice();
	mKneeCenter = femurKneeCenter;
	getNeckAxisVector();
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

Point HipFemur::getNeckCenter() const
{
	return mNeckCenter;
}

//Point HipFemur::getGreaterTrochanter() const
//{
//	return mGreaterTrochanter;
//}

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

Point HipFemur::getCanalAxisVectorInfSup() const
{
	return mCanalAxisVectorInfSup;
}

Point HipFemur::getVectorLatMed() const
{
	Line canalAxis(mCanalAxisVectorInfSup, mCanalAxisPoint);
	Point vector = mHeadCenter - canalAxis.getProjectPoint(mHeadCenter);
	vector.normalice();
	return vector;
}

Point HipFemur::getCanalAxisPoint() const
{
	return mCanalAxisPoint;
}

Point HipFemur::getNeckAxisVectorToHead() const
{
	return mNeckAxisVectorToHead;
}

Point HipFemur::getKneeCenter() const
{
	return mKneeCenter;
}

void HipFemur::getNeckAxisVector()
{
	/*Point midPoint = (mLesserTrochanter + mGreaterTrochanter) / 2.;

	Line canalAxis(mCanalAxisVectorInfSup, mCanalAxisPoint);

	midPoint = canalAxis.getProjectPoint(midPoint);

	Point tempAxis = mHeadCenter - midPoint;
	tempAxis.normalice();

	double radius = -1;
	Point neckCenter = midPoint;

	for (float i = 0.2; i < 0.9; i += 0.1)
	{
		Point cutPoint = midPoint + i * (mHeadCenter - midPoint);

		auto contour = ImplantTools::getContours(mFemur, tempAxis, cutPoint);
		if (contour->GetNumberOfPoints() != 0)
		{
			std::pair<Point, double> circle = ImplantTools::minCircle(contour, tempAxis);

			if (circle.second > 0 && (circle.second < radius || radius < 0))
			{
				radius = circle.second;
				neckCenter = circle.first;
			}
		}

	}*/

	mNeckAxisVectorToHead = mHeadCenter - mNeckCenter;
	mNeckAxisVectorToHead.normalice();
}

double HipFemur::getFemurVersion(const Point& pNeckAxisVectorToHead, const PelvisSide& pOperationSide) const
{
	Point neckAxisParameter = pNeckAxisVectorToHead;
	neckAxisParameter.normalice();

	Point lateralMedialAxis = mMedialEpicondyle - mLateralEpicondyle;
	lateralMedialAxis.normalice();
	Point midPoint = (mMedialEpicondyle + mLateralEpicondyle) / 2.;

	Plane axial;
	axial.init(mCanalAxisVectorInfSup, mCanalAxisPoint);

	Point neckAxis = axial.getProjectionVector(neckAxisParameter);
	lateralMedialAxis = axial.getProjectionVector(lateralMedialAxis);

	neckAxis.normalice();
	lateralMedialAxis.normalice();

	double angle = ImplantTools::getAngleBetweenVectors(neckAxis, lateralMedialAxis);

	Point anteriorDirectionVector;

	if (pOperationSide == PelvisSide::RIGHT_SIDE)
	{
		anteriorDirectionVector = lateralMedialAxis.cross(mCanalAxisVectorInfSup);
	}
	else
	{
		anteriorDirectionVector = mCanalAxisVectorInfSup.cross(lateralMedialAxis);
	}

	anteriorDirectionVector.normalice();

	Plane coronal;
	coronal.init(anteriorDirectionVector, midPoint);

	Point checkPoint = midPoint + 1000. * neckAxis;

	if (coronal.eval(checkPoint) >= 0)
	{
		return angle;
	}
	else
	{
		return -angle;
	}
}
