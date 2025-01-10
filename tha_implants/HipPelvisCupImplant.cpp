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

void HipPelvisCupImplant::init(const Point& pTopPoint, const Point& pBasePoint1, const Point& pBasePoint2, const Point& pBasePoint3, const std::vector<Point>& pExternalPoints, double pHemiSphereResolution, double pThickness)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_PELVIS_CUT_IMPLANT;
    }

    this->mTopPoint = pTopPoint;
    this->mBasePoint1 = pBasePoint1;
    this->mBasePoint2 = pBasePoint2;
    this->mBasePoint3 = pBasePoint3;
	this->mThickness = pThickness;

	mBasePlane.init(pBasePoint1, pBasePoint2, pBasePoint3);
	mBasePlane.reverseByPoint(pTopPoint);

	if (pExternalPoints.size() > 0)
	{
		auto fitSphere = ImplantTools::fitSphere(pExternalPoints);
		mCenter = fitSphere.first;
		mRadius = fitSphere.second;

		if (mRadius <= 0)
		{
			throw ImplantExceptionCode::CAN_NOT_FIT_SPHERE_TO_CUP_POINTS;
		}

		double distance = mBasePlane.getDistanceFromPoint(mCenter);
		double h = mRadius;
		if (distance < mRadius && distance > 0)
		{
			Point planeCutAndCenterVector = mCenter - mBasePlane.getProjectionPoint(mCenter);
			planeCutAndCenterVector.normalice();

			double direction = planeCutAndCenterVector.dot(mBasePlane.getNormalVector());
			
			if (direction > 0)
			{
				h = distance + mRadius;
			}
			else
			{
				h = mRadius - distance;
			}
		}

		mHemiSphereSurfaceArea = 2. * PI * mRadius * h;

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
		mRadius = -1;
	}

    isInit = true;
}

void HipPelvisCupImplant::init(const Point& pTopPoint, const Point& pBasePoint1, const Point& pBasePoint2, const Point& pBasePoint3, const Point& pSphereCenter, double pRadius, double pHemiSphereResolution, double pThickness)
{
	if (isInit == true)
	{
		throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_PELVIS_CUT_IMPLANT;
	}

	this->mTopPoint = pTopPoint;
	this->mBasePoint1 = pBasePoint1;
	this->mBasePoint2 = pBasePoint2;
	this->mBasePoint3 = pBasePoint3;
	this->mThickness = pThickness;

	mBasePlane.init(pBasePoint1, pBasePoint2, pBasePoint3);
	mBasePlane.reverseByPoint(pTopPoint);

	if (pRadius > 0)
	{
		mCenter = pSphereCenter;
		mRadius = pRadius;

		double distance = mBasePlane.getDistanceFromPoint(mCenter);
		double h = mRadius;
		if (distance < mRadius && distance > 0)
		{
			Point planeCutAndCenterVector = mCenter - mBasePlane.getProjectionPoint(mCenter);
			planeCutAndCenterVector.normalice();

			double direction = planeCutAndCenterVector.dot(mBasePlane.getNormalVector());

			if (direction > 0)
			{
				h = distance + mRadius;
			}
			else
			{
				h = mRadius - distance;
			}
		}

		mHemiSphereSurfaceArea = 2. * PI * mRadius * h;

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
		mRadius = -1;
	}

	isInit = true;
}

double HipPelvisCupImplant::getThickness() const
{
	return mThickness;
}

double HipPelvisCupImplant::getInternalRadius() const
{
	double result = mRadius - mThickness;
	if (result > 0)
	{
		return result;
	}
	return 0;
}

double HipPelvisCupImplant::getCupRadius() const
{
	return mRadius;
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
	if (mRadius <= 0)
	{
		Plane myPlane;
		myPlane.init(mBasePoint1, mBasePoint2, mBasePoint3);
		return myPlane.getProjectionPoint(mTopPoint);
	}
	else
	{
		return mCenter;
	}
}

Plane HipPelvisCupImplant::getBasePlane() const
{
	return mBasePlane;
}

vtkSmartPointer<vtkPolyData> HipPelvisCupImplant::getHemiSphereCup() const
{
	return mHemiSphereCup;
}

double HipPelvisCupImplant::getHemiSphereSurfaceArea() const
{
	return mHemiSphereSurfaceArea;
}