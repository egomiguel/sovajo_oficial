#include <iostream>
#include <math.h>
#include "RLine.hpp"
#include "RegistrationException.hpp"

using namespace UKA::REGISTRATION;

const double EPSILON = std::numeric_limits<double>::epsilon();

RLine::RLine(const cv::Point3d& directVector, const cv::Point3d& pPoint)
{
	this->mPoint = pPoint;
	this->directVector = directVector;
    normaliceDirectVector();
}

cv::Point3d RLine::getPoint() const { return mPoint; }

cv::Point3d RLine::getDirectVector() const { return directVector; }

double RLine::getSquareNorm(const cv::Point3d& pPoint) const
{
	return pPoint.x * pPoint.x + pPoint.y * pPoint.y + pPoint.z * pPoint.z;
}

double RLine::getDistanceFromPoint(const cv::Point3d& pPoint) const
{
	cv::Point3d hypotenuse = pPoint - mPoint;
	double tProjection = directVector.dot(hypotenuse) / getSquareNorm(directVector);
	cv::Point3d projection = directVector * tProjection;
	cv::Point3d distanceVector = hypotenuse - projection;
	return sqrt(getSquareNorm(distanceVector));
}

double RLine::getDistanceFromPoint(const double pPoint[3]) const
{
    cv::Point3d tPoint = cv::Point3d(pPoint[0], pPoint[1], pPoint[2]);
    return getDistanceFromPoint(tPoint);
}

RLine RLine::getParalellLine(const cv::Point3d& pPoint) const
{
	return RLine(directVector, pPoint);
}

RLine RLine::getPerpendicularLine(const cv::Point3d& pPoint) const
{
	if (isPointBelongToLine(pPoint) == true)
	{
        throw RegistrationExceptionCode::CAN_NOT_DETERMINE_PERPENDICULAR_LINE_TO_LINE_ON_REGISTRATION;
	}

	cv::Point3d hypotenuse = pPoint - mPoint;
	double tProjection = directVector.dot(hypotenuse) / getSquareNorm(directVector);
	cv::Point3d projection = directVector * tProjection;
	cv::Point3d perpendicularVector = hypotenuse - projection;

	return RLine(perpendicularVector, pPoint);
}

double RLine::getDistanceBetweenPoints(const cv::Point3d& a, const cv::Point3d& b, bool square)
{
	cv::Point3d diff = a - b;
	if (square == false)
	{
		return sqrt(diff.dot(diff));
	}
	else
	{
		return diff.dot(diff);
	}
}


cv::Point3d RLine::getPointAtDistance(const cv::Point3d& pPoint, const cv::Point3d& nearReferencePoint, float distance, bool closest) const
{
	cv::Point3d diffPoint = mPoint - pPoint;
	cv::Point3d equation;
	equation.x = directVector.dot(directVector);
	equation.y = 2.0*(directVector.dot(diffPoint));
	equation.z = (diffPoint.dot(diffPoint)) - (distance*distance);

	double disd = sqrt(equation.y*equation.y - 4.0 * equation.x * equation.z);
	double parameter1 = -(equation.y + disd) / (2.0 * equation.x);
	double parameter2 = -(equation.y - disd) / (2.0 * equation.x);
	cv::Point3d newPoint1 = (directVector * parameter1) + mPoint;
	cv::Point3d newPoint2 = (directVector * parameter2) + mPoint;
	cv::Point3d diff1 = nearReferencePoint - newPoint1;
	cv::Point3d diff2 = nearReferencePoint - newPoint2;
	if (closest == ((diff1.x * diff1.x + diff1.y * diff1.y + diff1.z * diff1.z) < (diff2.x * diff2.x + diff2.y * diff2.y + diff2.z * diff2.z)))
		return newPoint1;
	else
		return newPoint2;
}

cv::Point3d RLine::getProjectPoint(const cv::Point3d& pPoint) const
{
	cv::Point3d diff = pPoint - mPoint;
	cv::Point3d projectOnDirector = ((diff.dot(directVector)) / (directVector.dot(directVector))) * directVector;
	cv::Point3d projection = mPoint + projectOnDirector;
	return projection;
}

bool RLine::isPointBelongToLine(const cv::Point3d& pPoint) const
{
	double distance = getDistanceFromPoint(pPoint);
	if (distance < EPSILON)
	{
		return true;
	}
	return false;
}

cv::Mat RLine::getRotateMatrix(const cv::Point3d& axis, double angle)
{
    cv::Mat rotationMatrix(3, 3, CV_64F);

    cv::Point3d normaliceAxis = axis;
    normaliceAxis = normaliceAxis / sqrt(normaliceAxis.dot(normaliceAxis));
    
    rotationMatrix.at <double>(0, 0) = 0;
    rotationMatrix.at <double>(1, 0) = normaliceAxis.z;
    rotationMatrix.at <double>(2, 0) = -normaliceAxis.y;

    rotationMatrix.at <double>(0, 1) = -normaliceAxis.z;
    rotationMatrix.at <double>(1, 1) = 0;
    rotationMatrix.at <double>(2, 1) = normaliceAxis.x;

    rotationMatrix.at <double>(0, 2) = normaliceAxis.y;
    rotationMatrix.at <double>(1, 2) = -normaliceAxis.x;
    rotationMatrix.at <double>(2, 2) = 0;

    cv::Mat I = cv::Mat::eye(3, 3, CV_64F);
    cv::Mat result = I + sin(angle)*rotationMatrix + (1.0 - cos(angle))*(rotationMatrix * rotationMatrix);
    return result;
}

double RLine::getAngleBetweenVectors(const cv::Point3d& a, const cv::Point3d& b)
{
    double scalar = a.dot(b);
    double magnitude = sqrt((a.dot(a)) * (b.dot(b)));
    double tCos = scalar / magnitude;
    if (tCos <= -1.0)
    {
        return acos(-1.0);
    }
    else if (tCos >= 1.0)
    {
        return 0.0;
    }
    else
    {
        return acos(tCos);
    }
}

RLine RLine::makeLineWithPoints(const cv::Point3d& P1, const cv::Point3d& P2)
{
    cv::Point3d tempVector = P2 - P1;
    return RLine(tempVector, P1);
}

void RLine::normaliceDirectVector()
{
    directVector = directVector / sqrt(directVector.dot(directVector));
}

bool RLine::getInterceptionSphere(const cv::Point3d& center, double radius, std::pair<cv::Point3d, cv::Point3d>& result) const
{
    //ax^2 + bx + c = 0

    cv::Point3d constParam = mPoint - center;

    double a = directVector.dot(directVector);
    double b = 2.0 * (constParam.dot(directVector));
    double c = constParam.dot(constParam) - (radius * radius);

    double discriminant = b * b - 4.0 * a * c;

    if (discriminant < 0)
    {
        return false;
    }

    double x1 = (-b + sqrt(discriminant)) / (2.0 * a);
    double x2 = (-b - sqrt(discriminant)) / (2.0 * a);

    result.first = evalParameter(x1);
    result.second = evalParameter(x2);

    return true;
}

cv::Point3d RLine::evalParameter(double t) const
{
    double x = mPoint.x + (directVector.x * t);
    double y = mPoint.y + (directVector.y * t);
    double z = mPoint.z + (directVector.z * t);
    return cv::Point3d(x, y, z);
}