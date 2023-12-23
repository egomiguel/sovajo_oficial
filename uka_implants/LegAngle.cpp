#include <iostream>
#include "Utils.hpp"
#include "LegAngle.hpp"

using namespace UKA::IMPLANTS;

LegAngle::LegAngle()
{
	//isInit = false;
}


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
}