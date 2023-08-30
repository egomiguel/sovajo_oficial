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

void HipPelvisCupImplant::init(const Point& pTopPoint, const Point& pBasePoint1, const Point& pBasePoint2, const Point& pBasePoint3, const std::vector<Point>& pExternalPoints)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_PELVIS_CUT_IMPLANT;
    }

    this->mTopPoint = pTopPoint;
    this->mBasePoint1 = pBasePoint1;
    this->mBasePoint2 = pBasePoint2;
    this->mBasePoint3 = pBasePoint3;

	auto fitSphere = ImplantTools::fitSphere(pExternalPoints);
	Point center = fitSphere.first;
	double radius = fitSphere.second;
	double proportion = 0.5;

	if (radius == 0)
	{
		throw ImplantExceptionCode::CAN_NOT_FIT_SPHERE_TO_CUP_POINTS;
	}

	Plane basePlane;
	basePlane.init(pBasePoint1, pBasePoint2, pBasePoint3);
	basePlane.reverseByPoint(pTopPoint);

	double distance = basePlane.getDistanceFromPoint(center);
	if (distance < radius)
	{
		Point distal = center + radius * basePlane.getNormalVector();
		double distance = basePlane.getDistanceFromPoint(distal);
		proportion = distance / (2. * radius);
	}

	mHemiSphereSurfaceArea = proportion * 4. * PI * radius * radius;

	double pnt[3];
	pnt[0] = center.x;
	pnt[1] = center.y;
	pnt[2] = center.z;

	vtkNew<vtkSphereSource> sphere;
	sphere->SetCenter(pnt);
	sphere->SetRadius(radius);
	sphere->Update();

	vtkNew<vtkPlane> vtkPlaneA;

	auto planeNormal = basePlane.getNormalVector();
	auto planePoint = basePlane.getPoint();
	vtkPlaneA->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
	vtkPlaneA->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);

	vtkNew<vtkPlaneCollection> cutPlanes;
	cutPlanes->AddItem(vtkPlaneA);

	vtkNew<vtkClipClosedSurface> Clipper;
	Clipper->SetInputData(sphere->GetOutput());
	Clipper->SetClippingPlanes(cutPlanes);
	Clipper->Update();

	mHemiSphereCup = Clipper->GetOutput();

    isInit = true;
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