#include "HipPelvisCupImplant.hpp"
#include "ImplantsException.hpp"
#include "ImplantTools.hpp"
#include "vtkSphereSource.h"
#include "vtkPlane.h"
#include "vtkPlaneCollection.h"
#include "vtkClipClosedSurface.h"

using namespace THA::IMPLANTS;

HipPelvisCupImplant::HipPelvisCupImplant()
{
    isInit = false;
}

void HipPelvisCupImplant::init(const Point& pTopPoint, const Point& pBasePoint1, const Point& pBasePoint2, const Point& pBasePoint3, const std::vector<Point>& pExternalPoints, double pHemiSphereResolution)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_PELVIS_CUT_IMPLANT;
    }

    this->mTopPoint = pTopPoint;
    this->mBasePoint1 = pBasePoint1;
    this->mBasePoint2 = pBasePoint2;
    this->mBasePoint3 = pBasePoint3;

	if (pExternalPoints.size() > 0)
	{
		auto fitSphere = ImplantTools::fitSphere(pExternalPoints);
		mCenter = fitSphere.first;
		mRadius = fitSphere.second;
		double proportion = 0.5;

		if (mRadius <= 0)
		{
			throw ImplantExceptionCode::CAN_NOT_FIT_SPHERE_TO_CUP_POINTS;
		}

		mBasePlane.init(pBasePoint1, pBasePoint2, pBasePoint3);
		mBasePlane.reverseByPoint(pTopPoint);

		double distance = mBasePlane.getDistanceFromPoint(mCenter);
		if (distance < mRadius)
		{
			Point distal = mCenter + mRadius * mBasePlane.getNormalVector();
			double distance = mBasePlane.getDistanceFromPoint(distal);
			proportion = distance / (2. * mRadius);
		}

		mHemiSphereSurfaceArea = proportion * 4. * PI * mRadius * mRadius;

		double pnt[3];
		pnt[0] = mCenter.x;
		pnt[1] = mCenter.y;
		pnt[2] = mCenter.z;

		vtkNew<vtkSphereSource> sphere;

		double tResolution;
		double minResolution = sphere->GetThetaResolutionMinValue() / sphere->GetThetaResolutionMaxValue();

		if (pHemiSphereResolution <= 1 && pHemiSphereResolution >= minResolution)
		{
			tResolution = pHemiSphereResolution;
		}
		else if (pHemiSphereResolution > 1)
		{
			tResolution = 1;
		}
		else
		{
			tResolution = minResolution;
		}
	
		sphere->SetCenter(pnt);
		sphere->SetRadius(mRadius);
		sphere->SetThetaResolution(sphere->GetThetaResolutionMaxValue() * tResolution);
		sphere->SetPhiResolution(sphere->GetPhiResolutionMaxValue() * tResolution);
		sphere->Update();

		vtkNew<vtkPlane> vtkPlaneA;

		auto planeNormal = mBasePlane.getNormalVector();
		auto planePoint = mBasePlane.getPoint();
		vtkPlaneA->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
		vtkPlaneA->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);

		vtkNew<vtkPlaneCollection> cutPlanes;
		cutPlanes->AddItem(vtkPlaneA);

		vtkNew<vtkClipClosedSurface> Clipper;
		Clipper->SetInputData(sphere->GetOutput());
		Clipper->SetClippingPlanes(cutPlanes);
		Clipper->Update();

		mHemiSphereCup = Clipper->GetOutput();
	}
	else
	{
		mHemiSphereSurfaceArea = 0;
	}

    isInit = true;
}

void HipPelvisCupImplant::setHemisphereResolution(double pResolution)
{
	double pnt[3];
	pnt[0] = mCenter.x;
	pnt[1] = mCenter.y;
	pnt[2] = mCenter.z;

	vtkNew<vtkSphereSource> sphere;

	double tResolution;
	double minResolution = sphere->GetThetaResolutionMinValue() / sphere->GetThetaResolutionMaxValue();

	if (pResolution <= 1 && pResolution >= minResolution)
	{
		tResolution = pResolution;
	}
	else if (pResolution > 1)
	{
		tResolution = 1;
	}
	else
	{
		tResolution = minResolution;
	}

	sphere->SetCenter(pnt);
	sphere->SetRadius(mRadius);
	sphere->SetThetaResolution(sphere->GetThetaResolutionMaxValue() * tResolution);
	sphere->SetPhiResolution(sphere->GetPhiResolutionMaxValue() * tResolution);
	sphere->Update();

	vtkNew<vtkPlane> vtkPlaneA;

	auto planeNormal = mBasePlane.getNormalVector();
	auto planePoint = mBasePlane.getPoint();
	vtkPlaneA->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
	vtkPlaneA->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);

	vtkNew<vtkPlaneCollection> cutPlanes;
	cutPlanes->AddItem(vtkPlaneA);

	vtkNew<vtkClipClosedSurface> Clipper;
	Clipper->SetInputData(sphere->GetOutput());
	Clipper->SetClippingPlanes(cutPlanes);
	Clipper->Update();

	mHemiSphereCup = Clipper->GetOutput();
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

vtkSmartPointer<vtkPolyData> HipPelvisCupImplant::getHemiSphereCup() const
{
	return mHemiSphereCup;
}

double HipPelvisCupImplant::getHemiSphereSurfaceArea() const
{
	return mHemiSphereSurfaceArea;
}