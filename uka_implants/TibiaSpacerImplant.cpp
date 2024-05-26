#include <iostream>
#include "TibiaSpacerImplant.hpp"
#include "ImplantsException.hpp"

using namespace UKA::IMPLANTS;

TibiaSpacerImplant::TibiaSpacerImplant()
{
	isInit = false;
}

void TibiaSpacerImplant::init(const Point& apLinePclPoint, const Point& apLineTuberPoint, const Point& plateauRefPointUp, const Point& exteriorPointDown)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_TIBIA_SPACER_IMPLANT;
    }

	spacerPlane.init(apLinePclPoint, apLineTuberPoint, exteriorPointDown);
	kneeCenterPoint = (apLinePclPoint + apLineTuberPoint) / 2;

	Point newNormal = plateauRefPointUp - spacerPlane.getProjectionPoint(plateauRefPointUp);
	newNormal.normalice();
	spacerPlane.fixNormalVector(newNormal);

    this->plateauRefPointUp = plateauRefPointUp;
	this->pclPoint = apLinePclPoint;
	this->tuberPoint = apLineTuberPoint;
	this->plateauRefPointDown = exteriorPointDown;
    
    isInit = true;
}

TibiaSpacerImplant::TibiaSpacerImplant(const TibiaSpacerImplant& pImplant)
{
    this->spacerPlane = pImplant.spacerPlane;
    this->kneeCenterPoint = pImplant.kneeCenterPoint;
    this->isInit = pImplant.isInit;
	this->pclPoint = pImplant.pclPoint;
	this->tuberPoint = pImplant.tuberPoint;
	this->plateauRefPointDown = pImplant.plateauRefPointDown;
	this->plateauRefPointUp = pImplant.plateauRefPointUp;
}


Point TibiaSpacerImplant::getSpacerKneeCenter() const
{
    return kneeCenterPoint;
}

Plane TibiaSpacerImplant::getSpacerPlane() const
{
	return spacerPlane;
}

Point TibiaSpacerImplant::getSpacerNormalVector() const
{
	return spacerPlane.getNormalVector();
}

Point TibiaSpacerImplant::getSpacerVectorAP() const
{
	Point directVector = tuberPoint - pclPoint;
	directVector.normalice();
	return directVector;
}

Point TibiaSpacerImplant::getSpacerVectorTEA() const
{
    Point result = getSpacerNormalVector().cross(getSpacerVectorAP());
    result.normalice();
    return result;
}
