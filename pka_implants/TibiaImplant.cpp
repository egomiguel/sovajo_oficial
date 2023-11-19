
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
#include <iostream>
#include "TibiaImplant.hpp"
#include "ImplantsException.hpp"

using namespace PKA::IMPLANTS;

TibiaImplant::TibiaImplant()
{
	isInit = false;
}

void TibiaImplant::init(const Point& apLinePclPoint, const Point& apLineTuberPoint, const Point& sidePoint, const Point& exteriorPoint, const TibiaImplantInfo& pImplantInfo)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_TIBIA_IMPLANT;
    }
    this->mImplantInfo = pImplantInfo;
    tibiaPlane.init(apLinePclPoint, apLineTuberPoint, sidePoint);
    this->exteriorPoint = exteriorPoint;
	this->pclPoint = apLinePclPoint;
	this->tuberPoint = apLineTuberPoint;
	this->sidePoint = sidePoint;

    centralPoint = (apLinePclPoint + apLineTuberPoint) / 2;

    tibiaPlane.setPoint(centralPoint);
    Point newNormal = exteriorPoint - tibiaPlane.getProjectionPoint(exteriorPoint);
    tibiaPlane.fixNormalVector(newNormal);
    isInit = true;
}

TibiaImplant::TibiaImplant(const TibiaImplant& pImplant)
{
    this->tibiaPlane = pImplant.tibiaPlane;
    this->centralPoint = pImplant.centralPoint;
    this->isInit = pImplant.isInit;
    this->mImplantInfo = pImplant.mImplantInfo;
    this->exteriorPoint = pImplant.exteriorPoint;
	this->pclPoint = pImplant.pclPoint;
	this->tuberPoint = pImplant.tuberPoint;
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
	Point directVector = tuberPoint - pclPoint;
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

Point TibiaImplant::getPointPCL() const
{
	return pclPoint;
}

Point TibiaImplant::getPointTuber() const
{
	return tuberPoint;
}

Point TibiaImplant::getPointSide() const
{
	return sidePoint;
}
