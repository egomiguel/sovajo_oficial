
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
#include <iostream>
#include "TibiaImplant.hpp"
#include "ImplantsException.hpp"

using namespace TKA::IMPLANTS;

TibiaImplant::TibiaImplant()
{
	isInit = false;
}

//void TibiaImplant::fixNormalVectorTibia(const Point& fixPoint, const Point& referencePoint)
//{
//	Point tNormalVector = tibiaPlane.getNormalVector();
//	Point normalVector = Line::getFixDirectVector(tNormalVector, fixPoint, referencePoint, false);
//	tibiaPlane.fixNormalVector(normalVector);
//}

void TibiaImplant::init(const Point& pclPoint1, const Point& pclPoint2, const Point& frontPoint, const Point& exteriorPoint, const TibiaImplantInfo& pImplantInfo)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_TIBIA_IMPLANT;
    }
    this->mImplantInfo = pImplantInfo;
    tibiaPlane.init(pclPoint1, pclPoint2, frontPoint);
    this->exteriorPoint = exteriorPoint;

    midPoint = (pclPoint1 + pclPoint2) / 2;
    centralPoint = (frontPoint + midPoint) / 2;

    tibiaPlane.setPoint(centralPoint);
    Point newNormal = exteriorPoint - tibiaPlane.getProjectionPoint(exteriorPoint);
    tibiaPlane.fixNormalVector(newNormal);
    isInit = true;
}

TibiaImplant::TibiaImplant(const TibiaImplant& pImplant)
{
    this->tibiaPlane = pImplant.tibiaPlane;
    this->centralPoint = pImplant.centralPoint;
    this->midPoint = pImplant.midPoint;
    this->isInit = pImplant.isInit;
    this->mImplantInfo = pImplant.mImplantInfo;
    this->exteriorPoint = pImplant.exteriorPoint;
}

Point TibiaImplant::getExteriorPoint() const
{
    return exteriorPoint;
}

Point TibiaImplant::getCentralPoint() const
{
    return centralPoint;
}

Plane TibiaImplant::getTibiaPlane() const
{
	return tibiaPlane;
}

Point TibiaImplant::getTibiaNormalVector() const
{
	Point normalVector = tibiaPlane.getNormalVector();
	double squareNorm = normalVector.dot(normalVector);
	return normalVector / sqrt(squareNorm);
}

Point TibiaImplant::getTibiaVectorAP() const
{
	Point directVector = centralPoint - midPoint;
	double squareNorm = directVector.x * directVector.x + directVector.y * directVector.y + directVector.z * directVector.z;
	return directVector / sqrt(squareNorm);
}

Point TibiaImplant::getTibiaVectorTEA() const
{
    Point result = getTibiaNormalVector().cross(getTibiaVectorAP());
    result.normalice();
    return result;
}

cv::Mat TibiaImplant::getTibiaKneeCenter() const
{
	cv::Mat mat(3, 1, CV_64FC1);
	mat.at <double>(0, 0) = centralPoint.x;
	mat.at <double>(1, 0) = centralPoint.y;
	mat.at <double>(2, 0) = centralPoint.z;
	return mat;
}

TibiaImplantInfo TibiaImplant::getImplantInfo() const
{
    return mImplantInfo;
}

