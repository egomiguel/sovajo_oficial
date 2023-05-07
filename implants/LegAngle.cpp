#include <iostream>
#include "Utils.hpp"
#include "LegAngle.hpp"

using namespace TKA::IMPLANTS;

LegAngle::LegAngle()
{
	//isInit = false;
}

/*
void LegAngle::init(Point& hipCenterCT, Point& lateralEpicondyleCT, Point& medialEpicondyleCT, Point& femurKneeCenterCT)
{
	if (isInit == true)
	{
		std::cout << "Leg angle class has already been initialized." << std::endl;
		exit(1);
	}
	
	cv::Mat rotation(3, 3, CV_64FC1);
	cv::Mat translation(3, 1, CV_64FC1);

	for (int i = 0; i < 12; i += 4)
	{
		rotation.at<double>(i/4, 0) = TransformMatrixToReal[i + 0];
		rotation.at<double>(i/4, 1) = TransformMatrixToReal[i + 1];
		rotation.at<double>(i/4, 2) = TransformMatrixToReal[i + 2];

		translation.at<double>(i/4, 0) = TransformMatrixToReal[i + 3];
	}
	
	Point vectorAxis = hipCenterCT - femurKneeCenterCT;
	Point vectorML = lateralEpicondyleCT - medialEpicondyleCT;
	Point vectorAP = vectorAxis.cross(vectorML);

	coronalCT.init(vectorAP, hipCenterCT);
	this->hipCenterCT = hipCenterCT;
	
	isInit = true;
}
*/

cv::Mat LegAngle::pointToMat(const Point& pPoint) const
{
	cv::Mat mat(3, 1, CV_64FC1);
	mat.at <double>(0, 0) = pPoint.x;
	mat.at <double>(1, 0) = pPoint.y;
	mat.at <double>(2, 0) = pPoint.z;
	return mat;
}

double LegAngle::getAngleBetweenVectors(const Point& a, const Point& b) const
{
	double scalar = a.dot(b);
	double magnitude = sqrt((a.dot(a)) * (b.dot(b)));
	double tCos = scalar / magnitude;
	double angle;
	if (tCos <= -1.0)
	{
		angle = PI;
	}
	else if (tCos >= 1.0)
	{
		angle = 0;
	}
	else
	{
		angle = acos(tCos);
	}
	return (angle * 180.0) / PI;
}

double LegAngle::getFlexionAngle(const Point& hipCenter, const Point& femurKneeCenter, const Point& lateralEpicondyle, const Point& medialEpicondyle, const Point& tibiaKneeCenter, const Point& ankleCenter, bool isRight) const
{
    Point femurAxis = hipCenter - femurKneeCenter;
    Point tibiaAxis = ankleCenter - tibiaKneeCenter;
    Point femurTEA = lateralEpicondyle - medialEpicondyle;

    Plane sagital;
    sagital.init(femurTEA, femurKneeCenter);
    sagital.reverseByPoint(lateralEpicondyle);

    Point projA = sagital.getProjectionVector(femurAxis);
    Point projB = sagital.getProjectionVector(tibiaAxis);
    Point cross = projA.cross(projB);

    double flexion = 180.0 - abs(getAngleBetweenVectors(projA, projB));
    double sign = 1.0;

    Point ref = femurKneeCenter + 10.0 * cross;

    if (isRight == true)
    {
        if (sagital.eval(ref) < 0)
        {
            sign = -sign;
        }
    }
    else
    {
        if (sagital.eval(ref) > 0)
        {
            sign = -sign;
        }
    }

    return sign * flexion;

	/*double sign = 1.0;
    Point transverseVector = lateralEpicondyle - medialEpicondyle;
    Plane sagital;
    sagital.init(transverseVector, femurKneeCenter);

	Point directVectorFemur = hipCenter - femurKneeCenter;
	Point directVectorTibia = ankleCenter - tibiaKneeCenter;

    directVectorFemur = sagital.getProjectionVector(directVectorFemur);
    directVectorTibia = sagital.getProjectionVector(directVectorTibia);

	double angle = getAngleBetweenVectors(directVectorFemur, directVectorTibia);
	Line tibiaAxis((tibiaKneeCenter - ankleCenter), tibiaKneeCenter);
	Line perpendicularTibia = tibiaAxis.getPerpendicularLine(tibiaTubercle);

	Point tibiaHelpVector;

	Point a = tibiaTubercle + perpendicularTibia.getDirectVector();
	Point b = tibiaTubercle - perpendicularTibia.getDirectVector();

	if (tibiaAxis.getDistanceFromPoint(a) > tibiaAxis.getDistanceFromPoint(b))
	{
		tibiaHelpVector = a - tibiaTubercle;
	}
	else
	{
		tibiaHelpVector = b - tibiaTubercle;
	}

    tibiaHelpVector = sagital.getProjectionVector(tibiaHelpVector);

	if (angle > 100)
	{
		double helpAngle = getAngleBetweenVectors(tibiaHelpVector, directVectorFemur);
		if (helpAngle < 90.0)
			sign = -1.0;
		else
			sign = 1.0;
	}

	return (180.0 - angle)*sign;*/
}

double LegAngle::getVarusAngle(const Point& femurKneeCenter, const Point& medialCondyle, const Point& femurVectorAP, const Point& femurVectorTEA, const Point& tibiaVectorTEA, const Point& tibiaVectorAP, bool isRight) const
{
    Plane coronal, axialPlane, sagitalPlane;
    coronal.init(femurVectorAP, femurKneeCenter);
    coronal.reverseByPoint(medialCondyle, false);

	Point femurAxis, tibiaAxis;
	if (isRight == true)
	{
		femurAxis = femurVectorTEA.cross(femurVectorAP);
		tibiaAxis = tibiaVectorAP.cross(tibiaVectorTEA);
	}
	else
	{
		femurAxis = femurVectorAP.cross(femurVectorTEA);
		tibiaAxis = tibiaVectorTEA.cross(tibiaVectorAP);
	}

	axialPlane.init(femurAxis, femurKneeCenter);
	sagitalPlane.init(femurVectorTEA, femurKneeCenter);

	double flexion = abs(getAngleBetweenVectors(sagitalPlane.getProjectionVector(tibiaAxis), sagitalPlane.getProjectionVector(femurAxis)));

	Point projA, projB, cross;
	double varus;
	double sign = 1.0;

	if (flexion > 45 && flexion < 135)
	{
		projA = axialPlane.getProjectionVector(femurVectorTEA);
		projB = axialPlane.getProjectionVector(tibiaVectorTEA);
		cross = projA.cross(projB);

		varus = abs(getAngleBetweenVectors(projA, projB));

		Point ref = femurKneeCenter + 10.0 * cross;

		if (isRight == true)
		{
			if (axialPlane.eval(ref) < 0)
			{
				sign = -sign;
			}
		}
		else
		{
			if (axialPlane.eval(ref) > 0)
			{
				sign = -sign;
			}
		}
	}
	else
	{
		projA = coronal.getProjectionVector(femurVectorTEA);
		projB = coronal.getProjectionVector(tibiaVectorTEA);
		cross = projA.cross(projB);

		varus = abs(getAngleBetweenVectors(projA, projB));

		Point ref = femurKneeCenter + 10.0 * cross;

		if (isRight == true)
		{
			if (coronal.eval(ref) > 0)
			{
				sign = -sign;
			}
		}
		else
		{
			if (coronal.eval(ref) < 0)
			{
				sign = -sign;
			}
		}
	}

    //Inclination lateral is positive.

    
    return sign * varus;


	/*double sign;
	Plane coronal, tranverse, sagital;
    Point femurAxis = hipCenter - femurKneeCenter;

    tranverse.init(femurAxis, femurKneeCenter);

    Point vectorTEA = lateralEpicondyle - medialEpicondyle;
    Point vectorAP = vectorTEA.cross(femurAxis);

    coronal.init(vectorAP, femurKneeCenter);
    sagital.init(vectorTEA, femurKneeCenter);

    Point tibiaAxis = ankleCenter - tibiaKneeCenter;

    Side inclinationSide = getInclinationSide(hipCenter, femurKneeCenter, lateralEpicondyle, medialEpicondyle, ankleCenter);

	if (inclinationSide == LateralSide)
	{
		sign = 1.0;
	}
	else
	{
		sign = -1.0;
	}

    double flexion = abs(getAngleBetweenVectors(sagital.getProjectionVector(femurAxis), sagital.getProjectionVector(tibiaAxis)));
    double angle;

    if ( flexion > 45.0 && flexion < 135.0 )
    {
        Point coronalNormalProj = tranverse.getProjectionVector(coronal.getNormalVector());
        Point tibiaAxisProj = tranverse.getProjectionVector(tibiaAxis);
        angle = abs(getAngleBetweenVectors(coronalNormalProj, tibiaAxisProj));
        if (angle > 90.0)
        {
            angle = 180.0 - angle;
        }
    }
    else
    {
        Point femurProjectionVector = coronal.getProjectionVector(femurAxis);
        Point tibiaProjectionVector = coronal.getProjectionVector(tibiaAxis);
        angle = abs(getAngleBetweenVectors(femurProjectionVector, tibiaProjectionVector));

        if (angle > 90.0)
        {
            angle = 180.0 - angle;
        }
    }
	return sign * angle;*/

}

double LegAngle::getRotationAngle(const Point& hipCenter, const Point& femurKneeCenter, const Point& femurVectorTEA, const Point& tibiaVectorTEA, const Point& tibiaVectorAP, const Point& unkle, bool isRight) const
{
    Plane axialPlane, axialPlaneTibia, sagitalPlane, coronalPlane;
    Point vectorAxis = hipCenter - femurKneeCenter;
	Point vectorAxisTibia = tibiaVectorTEA.cross(tibiaVectorAP);

    axialPlane.init(vectorAxis, femurKneeCenter);
    axialPlane.reverseByPoint(hipCenter);

	axialPlaneTibia.init(vectorAxisTibia, femurKneeCenter);
	axialPlaneTibia.reverseByPoint(unkle);
	vectorAxisTibia = axialPlaneTibia.getNormalVector();

	sagitalPlane.init(femurVectorTEA, femurKneeCenter);

	Point femurAP;
	if (isRight == true)
	{
		femurAP = axialPlane.getNormalVector().cross(femurVectorTEA);
	}
	else
	{
		femurAP = femurVectorTEA.cross(axialPlane.getNormalVector());
	}

	coronalPlane.init(femurAP, femurKneeCenter);

	double flexion = abs(getAngleBetweenVectors(sagitalPlane.getProjectionVector(vectorAxisTibia), sagitalPlane.getProjectionVector(axialPlane.getNormalVector())));

	Point projA, projB, cross;
	double rotation;
	double sign = 1.0;

	if (flexion > 45 && flexion < 135)
	{
		projA = coronalPlane.getProjectionVector(femurVectorTEA);
		projB = coronalPlane.getProjectionVector(tibiaVectorTEA);
		cross = projA.cross(projB);
		rotation = abs(getAngleBetweenVectors(projA, projB));

		Point ref = femurKneeCenter + 10.0 * cross;

		if (isRight == true)
		{
			if (coronalPlane.eval(ref) < 0)
			{
				sign = -sign;
			}
		}
		else
		{
			if (coronalPlane.eval(ref) > 0)
			{
				sign = -sign;
			}
		}
	}
	else
	{
		projA = axialPlane.getProjectionVector(femurVectorTEA);
		projB = axialPlane.getProjectionVector(tibiaVectorTEA);
		cross = projA.cross(projB);
		rotation = abs(getAngleBetweenVectors(projA, projB));

		Point ref = femurKneeCenter + 10.0 * cross;

		if (isRight == true)
		{
			if (axialPlane.eval(ref) < 0)
			{
				sign = -sign;
			}
		}
		else
		{
			if (axialPlane.eval(ref) > 0)
			{
				sign = -sign;
			}
		}
	}
    
    //If femur vector rotate medial to lateral respect tibia vector, then is positive.


    return sign * rotation;
}

Side LegAngle::getInclinationSide(const Point& hipCenter, const Point& kneeCenter, const Point& latEpi, const Point& medEpi, const Point& ankle) const
{
    Point femurAxis = hipCenter - kneeCenter;
    Point vectorTEA = latEpi - medEpi;
    //Point vectorAP = vectorTEA.cross(femurAxis);
    Plane sagital, transverse;
    transverse.init(femurAxis, kneeCenter);
    sagital.init(vectorTEA, kneeCenter);
    Point ankleProj = transverse.getProjectionPoint(ankle);

    if (sagital.eval(latEpi) < 0)
    {
        sagital.reverse();
    }

    double eval = sagital.eval(ankleProj);

    if (eval > 0)
    {
        return LateralSide;
    }
    else
    {
        return MedialSide;
    }

    /*Point femurAxis = hipCenter - kneeCenter;
    Point legAxis = hipCenter - ankle;
    Plane transverse, coronal;
    transverse.init(femurAxis, kneeCenter);

    Point vectorTEA = latEpi - medEpi;
    Point vectorAP = vectorTEA.cross(femurAxis);

    coronal.init(vectorAP, kneeCenter);

    Point hipProj = coronal.getProjectionPoint(hipCenter);
    Point latProj = coronal.getProjectionPoint(latEpi);
    Point medProj = coronal.getProjectionPoint(medEpi);
    Point kneeProj = coronal.getProjectionPoint(kneeCenter);
    Point legAxisProj = coronal.getProjectionVector(legAxis);

    Line lineRef = Line(legAxisProj, hipProj);
    Point epiVector = latProj - medProj;
    epiVector.normalice();

    Point latRef = kneeProj + epiVector;
    Point medRef = kneeProj - epiVector;

    if (lineRef.getDistanceFromPoint(latRef) > lineRef.getDistanceFromPoint(medRef))
    {
        return MedialSide;
    }
    else
    {
        return LateralSide;
    }*/
}