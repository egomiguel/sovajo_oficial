#include <iostream>
#include <math.h>
#include "Utils.hpp"
#include "Line.hpp"
#include "ImplantsException.hpp"

using namespace TKA::IMPLANTS;

Line::Line(const Point& directVector, const Point& pPoint)
{
    if (directVector.dot(directVector) == 0)
    {
        throw ImplantExceptionCode::LINE_DIRECT_VECTOR_CAN_NOT_BE_ZERO;
    }
	this->mPoint = pPoint;
	this->directVector = directVector;
    normaliceDirectVector();
}

Line::Line(const Line& pLine)
{
    this->directVector = pLine.getDirectVector();
    this->mPoint = pLine.getPoint();
}

void Line::setDirectVector(const Point& newVector)
{
    if (newVector.dot(newVector) == 0)
    {
        throw ImplantExceptionCode::LINE_DIRECT_VECTOR_CAN_NOT_BE_ZERO;
    }
    this->directVector = newVector;
    normaliceDirectVector();
}

Point Line::getPoint() const { return mPoint; }

Point Line::getDirectVector() const { return directVector; }

double Line::getSquareNorm(const Point& pPoint) const
{
	return pPoint.x * pPoint.x + pPoint.y * pPoint.y + pPoint.z * pPoint.z;
}

double Line::getDistanceFromPoint(const Point& pPoint) const
{
	Point hypotenuse = pPoint - mPoint;
	double tProjection = directVector.dot(hypotenuse) / getSquareNorm(directVector);
	Point projection = directVector * tProjection;
	Point distanceVector = hypotenuse - projection;
	return sqrt(getSquareNorm(distanceVector));
}

Line Line::getParalellLine(const Point& pPoint) const
{
	return Line(directVector, pPoint);
}

Line Line::getPerpendicularLine(const Point& pPoint) const
{
	if (isPointBelongToLine(pPoint) == true)
	{
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_PERPENDICULAR_LINE_TO_LINE;
	}

	Point hypotenuse = pPoint - mPoint;
	double tProjection = directVector.dot(hypotenuse) / getSquareNorm(directVector);
	Point projection = directVector * tProjection;
	Point perpendicularVector = hypotenuse - projection;

	return Line(perpendicularVector, pPoint);
}

bool Line::isPointBelongToLine(const Point& pPoint) const
{
	double distance = getDistanceFromPoint(pPoint);
	if (distance < EPSILON)
	{
		return true;
	}
	return false;
}

Point Line::getFixDirectVector(const Point& pLineDirectVector, const Point& pLineFixPoint, const Point& pReferencePoint, bool closest)
{
	Point tPoint1 = pLineFixPoint + pLineDirectVector;
	Point tPoint2 = pLineFixPoint - pLineDirectVector;
	Point diff1 = pReferencePoint - tPoint1;
	Point diff2 = pReferencePoint - tPoint2;

    if (closest == true)
    {
        if ((diff1.x * diff1.x + diff1.y * diff1.y + diff1.z * diff1.z) < (diff2.x * diff2.x + diff2.y * diff2.y + diff2.z * diff2.z))
            return (tPoint1 - pLineFixPoint);
        else
            return (tPoint2 - pLineFixPoint);
    }
    else
    {
        if ((diff1.x * diff1.x + diff1.y * diff1.y + diff1.z * diff1.z) < (diff2.x * diff2.x + diff2.y * diff2.y + diff2.z * diff2.z))  
            return (tPoint2 - pLineFixPoint);
        else
            return (tPoint1 - pLineFixPoint);
    }

}

cv::Mat Line::getRotateMatrix(const Point& axis, double angle)
{
	cv::Mat rotationMatrix(3, 3, CV_64F);

    Point normaliceAxis = axis;
    normaliceAxis.normalice();

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

double Line::getAngleBetweenVectors(const Point& a, const Point& b)
{
	double scalar = a.dot(b);
	double magnitude = sqrt((a.dot(a)) * (b.dot(b)));
	double tCos = scalar / magnitude;
	if (tCos <= -1.0)
	{
		return PI;
	}
	else if (tCos >= 1.0)
	{
		return 0;
	}
	else
	{
		return acos(tCos);
	}
}

double Line::getDistanceBetweenPoints(const Point& a, const Point& b, bool square)
{
	Point diff = a - b;
	if (square == false)
	{
		return sqrt(diff.dot(diff));
	}
	else
	{
		return diff.dot(diff);
	}
}

Point Line::getPointAtDistance(float distance) const
{
	return mPoint + distance * directVector;
}

Point Line::getPointAtDistance(const Point& pPoint, const Point& nearReferencePoint, float distance, bool closest) const
{
	/*Point diffPoint = mPoint - pPoint;
	Point equation;
	equation.x = directVector.dot(directVector);
	equation.y = 2.0*(directVector.dot(diffPoint));
	equation.z = (diffPoint.dot(diffPoint)) - (distance*distance);
	double disd = sqrt(equation.y*equation.y - 4.0 * equation.x * equation.z);
	double parameter1 = -(equation.y + disd) / (2.0 * equation.x);
	double parameter2 = -(equation.y - disd) / (2.0 * equation.x);
	Point newPoint1 = (directVector * parameter1) + mPoint;
	Point newPoint2 = (directVector * parameter2) + mPoint;*/

    Point myDirectVector = directVector;
    myDirectVector.normalice();

    Point newPoint1 = pPoint + (distance * myDirectVector);
    Point newPoint2 = pPoint - (distance * myDirectVector);

	Point diff1 = nearReferencePoint - newPoint1;
	Point diff2 = nearReferencePoint - newPoint2;
	if (closest == ((diff1.x * diff1.x + diff1.y * diff1.y + diff1.z * diff1.z) < (diff2.x * diff2.x + diff2.y * diff2.y + diff2.z * diff2.z)))
		return newPoint1;
	else
		return newPoint2;
}

Point Line::getProjectPoint(const Point& pPoint) const
{
	Point diff = pPoint - mPoint;
	Point projectOnDirector = ((diff.dot(directVector)) / (directVector.dot(directVector))) * directVector;
	Point projection = mPoint + projectOnDirector;
	return projection;
}

Point Line::getProjectPoint(const Point& linePoint1, const Point& linePoint2, const Point& externalPoint)
{
    Point tDirecVector = linePoint1 - linePoint2;
    Point diff = externalPoint - linePoint1;
    Point projectOnDirector = ((diff.dot(tDirecVector)) / (tDirecVector.dot(tDirecVector))) * tDirecVector;
    Point projection = linePoint1 + projectOnDirector;
    return projection;
}


void Line::setPoint(const Point& newPoint)
{
    this->mPoint = newPoint;
}

void Line::normaliceDirectVector()
{
    directVector = directVector / sqrt(directVector.dot(directVector));
}

Line Line::makeLineWithPoints(const Point& P1, const Point& P2)
{

    Point tempVector = P2 - P1;
    return Line(tempVector, P1);
}

Line Line::getBestLine(const std::vector<Point>& pPoints)
{
	if (pPoints.size() < 2)
	{
		throw ImplantExceptionCode::LINE_DIRECT_VECTOR_CAN_NOT_BE_ZERO;
	}

	cv::Point3d centroid(0, 0, 0);
	for (const auto& pt : pPoints) {
		centroid += pt;
	}
	centroid *= (1.0 / pPoints.size());

	cv::Matx33d covariance = cv::Matx33d::zeros();
	for (const auto& pt : pPoints) {
		cv::Vec3d diff = cv::Vec3d(pt.x - centroid.x, pt.y - centroid.y, pt.z - centroid.z);
		covariance += diff * diff.t();
	}

	cv::SVD svd(covariance);
	Point directVector = Point(svd.u.at<double>(0, 0), svd.u.at<double>(1, 0), svd.u.at<double>(2, 0));
	Point pointOnLine = centroid;

	return Line(directVector, pointOnLine);
}

bool Line::getInterceptionSphere(const Point& center, double radius, std::pair<Point, Point>& result) const
{
    //ax^2 + bx + c = 0

    Point constParam = mPoint - center;

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

Point Line::evalParameter(double t) const
{
    double x = mPoint.x + (directVector.x * t);
    double y = mPoint.y + (directVector.y * t);
    double z = mPoint.z + (directVector.z * t);
    return Point(x, y, z);
}

Point Line::getNearestPoint(const std::vector<Point>& points)
{
    auto it1 = points.begin();
    auto it2 = points.end();

    double dist = -1;
    Point nearPoint;

    for ( ; it1 != it2; ++it1)
    {
        double temp = getDistanceFromPoint(*it1);
        
        if (temp < dist || dist < 0)
        {
            dist = temp;
            nearPoint = *it1;
        }
    }
    return nearPoint;
}

