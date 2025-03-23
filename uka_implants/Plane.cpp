#include <iostream>
#include <math.h>
#include <iomanip>
#include "Utils.hpp"
#include "Plane.hpp"
#include "ImplantsException.hpp"

using namespace UKA::IMPLANTS;

struct PointsComp
{
	cv::Mat rotation;
    Point centroid;
	PointsComp(const cv::Mat& pRotation, Point pCentroid) : rotation(pRotation), centroid(pCentroid)
    {
        /*
        cv::Mat pointMat(3, 1, CV_64F);

        pointMat.at <double>(0, 0) = pCentroid.x;
        pointMat.at <double>(1, 0) = pCentroid.y;
        pointMat.at <double>(2, 0) = pCentroid.z;

        cv::Mat newCentroid = pRotation * pointMat;
        centroid = Point(newCentroid);*/
    }

	bool operator()(cv::Point3d pPoint1, cv::Point3d pPoint2) const
	{
		cv::Mat pointMat1(3, 1, CV_64F);
		cv::Mat pointMat2(3, 1, CV_64F);

        cv::Point3d point1 = pPoint1 - centroid;
        cv::Point3d point2 = pPoint2 - centroid;

		pointMat1.at <double>(0, 0) = point1.x;
		pointMat1.at <double>(1, 0) = point1.y;
		pointMat1.at <double>(2, 0) = point1.z;

		pointMat2.at <double>(0, 0) = point2.x;
		pointMat2.at <double>(1, 0) = point2.y;
		pointMat2.at <double>(2, 0) = point2.z;

		cv::Mat rotateP1 = rotation * pointMat1;
		cv::Mat rotateP2 = rotation * pointMat2;

		double angle_p1 = atan2(rotateP1.at<double>(1, 0), rotateP1.at<double>(0, 0)) * 180.0 / PI;
		double angle_p2 = atan2(rotateP2.at<double>(1, 0), rotateP2.at<double>(0, 0)) * 180.0 / PI;
		if (angle_p1 < 0)
		{
			angle_p1 = 360.0 + angle_p1;
		}
		if (angle_p2 < 0)
		{
			angle_p2 = 360.0 + angle_p2;
		}
		return (angle_p1 < angle_p2);
	}
};

Plane::Plane()
{
	normalVector = Point(0.0, 0.0, 0.0);
	bias = 0.0;
	isInit = false;
}

void Plane::fixNormalVector(const Point& newNormalVector)
{
	normalVector = newNormalVector;
	bias = (-1.0)*newNormalVector.dot(mPoint);
    normalizeNormalVector();
}

void Plane::setPoint(const Point& pPoint)
{
	mPoint = pPoint;
	bias = (-1.0)*normalVector.dot(pPoint);
}

void Plane::movePlane(const Point& pPoint)
{
    mPoint = pPoint;
    bias = (-1.0)*normalVector.dot(pPoint);
}

void Plane::init(const Point& point1, const Point& point2, const Point& point3)
{
    if (point1 == point2 || point1 == point3 || point2 == point3)
    {
        throw ImplantExceptionCode::POINTS_TO_DEFINE_PLANE_CAN_NOT_BE_EQUAL;
    }

	Line temp((point1 - point2), point2);
	if (temp.isPointBelongToLine(point3))
	{
        throw ImplantExceptionCode::POINTS_TO_DEFINE_PLANE_CAN_NOT_BE_COLLINEAR;
	}
	Point directVector1, directVector2;
	if (isInit == true)
	{
        throw ImplantExceptionCode::ALREADY_INITIALIZED_PLANE;
	}

	directVector1 = point2 - point1;
	directVector2 = point3 - point1;

	mPoint = point3;
	normalVector = directVector1.cross(directVector2);
	bias = (-1.0)*normalVector.dot(point3);
    normalizeNormalVector();
	isInit = true;
}

void Plane::init(const Point& normalVector, const Point& pPoint)
{
	if (isInit == true)
	{
        throw ImplantExceptionCode::ALREADY_INITIALIZED_PLANE;
	}

    if (normalVector.dot(normalVector) == 0)
    {
        throw ImplantExceptionCode::PLANE_NORMAL_VECTOR_CAN_NOT_BE_ZERO;
    }

	mPoint = pPoint;
	this->normalVector = normalVector;
	this->bias = (-1.0) * normalVector.dot(pPoint);
    normalizeNormalVector();
	isInit = true;
}

void Plane::deletePlane()
{
    normalVector = Point(0.0, 0.0, 0.0);
    mPoint = Point(0.0, 0.0, 0.0);
    bias = 0.0;
    isInit = false;
}

bool Plane::getIsInit() const
{
    return isInit;
}

Plane Plane::getBestPlane(const std::vector<Point>& pPoints, bool& result)
{
    Plane bestPlane;
    int tSize = pPoints.size();
    cv::Point3d averagePoint(0, 0, 0);
    std::vector<cv::Point3d> myPoints;

    for (int i = 0; i < tSize; i++)
    {
        cv::Point3d pnt = { pPoints[i].x, pPoints[i].y, pPoints[i].z };
        myPoints.push_back(pnt);
        averagePoint = averagePoint + cv::Point3d(pPoints[i].x, pPoints[i].y, pPoints[i].z);
    }

    if (tSize > 2)
    {
        averagePoint = averagePoint / double(tSize);
    }
    else
    {
        result = false;
        return bestPlane;
    }

    std::vector<cv::Point3d> center(tSize, averagePoint);

    cv::Mat A = cv::Mat(tSize, 3, CV_64F, myPoints.data());
    cv::Mat C(tSize, 3, CV_64F, center.data());

    //Subtracting centroid_ from the set of points to center on the origin.
    cv::Mat CA = A - C;
    cv::SVD svd(CA);
    cv::Mat Last_V_Row = svd.vt.row(2);
    double a = Last_V_Row.at<double>(0);
    double b = Last_V_Row.at<double>(1);
    double c = Last_V_Row.at<double>(2);

    bestPlane.init(Point(a, b, c), averagePoint);
    result = true;
    //bestPlane.show();
    return bestPlane;
}

Plane Plane::getBestPlane(const std::list<Point>& pPoints, bool& result)
{
    Plane bestPlane;
    int tSize = pPoints.size();
    cv::Point3d averagePoint(0, 0, 0);
    std::vector<cv::Point3d> myPoints;

    auto it1 = pPoints.begin();
    auto it2 = pPoints.end();

    for ( ; it1 != it2; ++it1)
    {
        myPoints.push_back((*it1).ToCVPoint());
        averagePoint = averagePoint + (*it1).ToCVPoint();
    }

    if (tSize > 2)
    {
        averagePoint = averagePoint / double(tSize);
    }
    else
    {
        result = false;
        return bestPlane;
    }

    std::vector<cv::Point3d> center(tSize, averagePoint);

    cv::Mat A = cv::Mat(tSize, 3, CV_64F, myPoints.data());
    cv::Mat C(tSize, 3, CV_64F, center.data());

    //Subtracting centroid_ from the set of points to center on the origin.
    cv::Mat CA = A - C;
    cv::SVD svd(CA);
    cv::Mat Last_V_Row = svd.vt.row(2);
    double a = Last_V_Row.at<double>(0);
    double b = Last_V_Row.at<double>(1);
    double c = Last_V_Row.at<double>(2);

    bestPlane.init(Point(a, b, c), averagePoint);
    result = true;
    //bestPlane.show();
    return bestPlane;
}

Plane Plane::getBestPlaneOutliers(const std::vector<Point>& pPoints, std::vector<Point>& pOutLiers, bool& result, double pPercent)
{
    Plane bestPlane;
    int tSize = pPoints.size();

    if (tSize < 3 || pPercent > 0.5)
    {
        result = false;
        return bestPlane;
    }
    else if (tSize == 3)
    {
        result = true;
        bestPlane.init(pPoints[0], pPoints[1], pPoints[2]);
        return bestPlane;
    }
    
    pOutLiers.clear();
    std::list<Point> inLiers;
    int outAmoung = tSize * pPercent;

    for (int i = 0; i < tSize; i++)
    {
        inLiers.push_back(pPoints[i]);
    }

    while (pOutLiers.size() <= outAmoung)
    {
        bestPlane = Plane::getBestPlane(inLiers, result);
        if (result == false)
        {
            break;
        }

        if (pOutLiers.size() == outAmoung)
        {
            break;
        }

        auto it1 = inLiers.begin();
        auto it2 = inLiers.end();

        double distance = -1;
        std::list<Point>::iterator it;

        for ( ; it1 != it2; ++it1)
        {
            double temp = abs(bestPlane.eval(*it1));
            if (temp > distance)
            {
                distance = temp;
                it = it1;
            }
        }

        pOutLiers.push_back(*it);
        inLiers.erase(it);
    }

    return bestPlane;
}

Plane Plane::getBestPlane(const std::vector<PointTypeITK>& pPoints, bool& result)
{
    Plane bestPlane;
    std::vector<cv::Point3d> myPoints;
    int tSize = pPoints.size();
    cv::Point3d averagePoint(0, 0, 0);

    for (int i = 0; i < tSize; i++)
    {
        cv::Point3d pnt = { pPoints[i][0], pPoints[i][1], pPoints[i][2] };
        myPoints.push_back(pnt);

        averagePoint = averagePoint + cv::Point3d(pPoints[i][0], pPoints[i][1], pPoints[i][2]);
    }

    if (tSize > 0)
    {
        averagePoint = averagePoint / double(tSize);
    }
    else
    {
        result = false;
        return bestPlane;
    }

    std::vector<cv::Point3d> center(tSize, averagePoint);

    cv::Mat A = cv::Mat(tSize, 3, CV_64F, myPoints.data());
    cv::Mat C(tSize, 3, CV_64F, center.data());

    //Subtracting centroid_ from the set of points to center on the origin.
    cv::Mat CA = A - C;
    cv::SVD svd(CA);
    cv::Mat Last_V_Row = svd.vt.row(2);
    double a = Last_V_Row.at<double>(0);
    double b = Last_V_Row.at<double>(1);
    double c = Last_V_Row.at<double>(2);

    bestPlane.init(Point(a, b, c), averagePoint);
    result = true;
    //bestPlane.show();
    return bestPlane;
}

Plane Plane::getBestPlane(const std::vector<Point>& pPoints, const Point& pParallelVector, bool& result)
{
    throw ("No implemented yet");

    Plane bestPlane;
    double x = pParallelVector.x;
    double y = pParallelVector.y;
    double z = pParallelVector.z;

    if (abs(x) == abs(y) && abs(y) == abs(z) && abs(x) == 0)
    {
        result = false;
        return bestPlane;
    }

    if (abs(x) >= abs(y) && abs(y) >= abs(z))
    {


    }
    else if (abs(y) >= abs(x) && abs(x) >= abs(z))
    {

    }
    else
    {

    }
}

double Plane::getSquareNorm(const Point& pPoint) const
{
	return pPoint.x * pPoint.x + pPoint.y * pPoint.y + pPoint.z * pPoint.z;
}

void Plane::normalizeNormalVector()
{
	double squareNorm = sqrt(normalVector.dot(normalVector));
	normalVector = (normalVector / squareNorm);
	bias = (bias / squareNorm);
}

void Plane::transformPlane(const cv::Mat& rotation, const cv::Mat& translation)
{
	cv::Mat transformNormalVector = rotation * normalVector.ToMatPoint();
	cv::Mat transformPoint = (rotation * mPoint.ToMatPoint()) + translation;

	mPoint = Point(transformPoint);
	this->normalVector = Point(transformNormalVector);
	this->bias = (-1.0) * normalVector.dot(mPoint);
	normalizeNormalVector();
}

cv::Mat Plane::getNormalVectorMat() const
{
	cv::Mat mat(3, 1, CV_64F);
	mat.at <double>(0, 0) = normalVector.x;
	mat.at <double>(1, 0) = normalVector.y;
	mat.at <double>(2, 0) = normalVector.z;
	return mat;
}

cv::Mat Plane::getPointMat() const
{
	cv::Mat mat(3, 1, CV_64FC1);
	mat.at <double>(0, 0) = mPoint.x;
	mat.at <double>(1, 0) = mPoint.y;
	mat.at <double>(2, 0) = mPoint.z;
	return mat;
}

Plane Plane::getPerpendicularPlane(const Point& point1, const Point& point2) const
{
	Point tDirectVector = point2 - point1;
	Point tNormalVector = normalVector.cross(tDirectVector);

	if (tNormalVector.dot(tNormalVector) < EPSILON)
	{
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_PERPENDICULAR_PLANE_TO_PLANE;
	}

	Plane perpendicularPlane;
	perpendicularPlane.init(tNormalVector, point1);
	return perpendicularPlane;
}

Point Plane::getInterceptionLinePoint(const Line& a) const
{
	double magnitude = normalVector.dot(a.getDirectVector());
	if (magnitude == 0.0)
	{
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_INTERCEPTION_LINE_TO_PLANE;
	}
	double parameter = -(bias + (normalVector.dot(a.getPoint()))) / magnitude;
	
	Point result = a.getPoint() + (a.getDirectVector() * parameter);
	
	return result;
}

Point Plane::getProjectionVector(const Point& vector) const
{
	Point projectOnNormal = ((vector.dot(normalVector)) / (normalVector.dot(normalVector))) * normalVector;
	Point projection = vector - projectOnNormal;
	return projection;
}

Point Plane::getProjectionPoint(const Point& pPoint) const
{
	Point diff = pPoint - mPoint;
	Point projectOnNormal = ((diff.dot(normalVector)) / (normalVector.dot(normalVector))) * normalVector;
	Point projection = pPoint - projectOnNormal;
	return projection;
}

bool Plane::isPointBelongToPlane(const Point& pPoint) const
{
	double evaluation = normalVector.dot(pPoint) + bias;
	if (fabs(evaluation) < EPSILON)
	{
		return true;
	}
	return false;
}

bool Plane::isPointNearToPlane(const Point& pPoint, double distance) const
{
	double evaluation = (normalVector.dot(pPoint) + bias) / sqrt(getSquareNorm(normalVector));
	if (fabs(evaluation) <= distance)
	{
		return true;
	}
	return false;
}

double Plane::getDistanceFromPoint(const Point& pPoint) const
{
	double distance = abs(normalVector.dot(pPoint) + bias) / sqrt(getSquareNorm(normalVector));
	return distance;
}

Point Plane::getNormalVector() const
{
	return normalVector;
}

double Plane::getBias() const
{
	return bias;
}

Point Plane::getPoint() const
{
	return mPoint;
}

double Plane::eval(const Point& a) const
{
	return normalVector.dot(a) + bias;
}

void Plane::sortCoplanarPointsByAngle(std::vector<Point>& points, Point& centroid, Point& fixPlaneNormalIn, bool clockwise)
{
	Point normalXY(0.0, 0.0, 1.0);
	fixPlaneNormalIn = fixPlaneNormalIn / sqrt(fixPlaneNormalIn.dot(fixPlaneNormalIn));
	Point rotationAxis = normalXY.cross(fixPlaneNormalIn);
	rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
	double rotationAngle = Line::getAngleBetweenVectors(normalXY, fixPlaneNormalIn);
	cv::Mat rotation_1 = Line::getRotateMatrix(rotationAxis, -rotationAngle);
	cv::Mat rotation_2 = Line::getRotateMatrix(rotationAxis, rotationAngle);

	cv::Mat normalMat(3, 1, CV_64F);
	normalMat.at<double>(0, 0) = fixPlaneNormalIn.x;
	normalMat.at<double>(1, 0) = fixPlaneNormalIn.y;
	normalMat.at<double>(2, 0) = fixPlaneNormalIn.z;

	cv::Mat rotateVector_1 = rotation_1 * normalMat;
	cv::Mat rotateVector_2 = rotation_2 * normalMat;

	cv::Mat rotate;

	double distance_1 = Line::getDistanceBetweenPoints(Point(rotateVector_1), normalXY);
	double distance_2 = Line::getDistanceBetweenPoints(Point(rotateVector_2), normalXY);

	if (distance_1 < distance_2)
	{
		rotate = rotation_1;
	}
	else
	{
		rotate = rotation_2;
	}
	sort(points.begin(), points.end(), PointsComp(rotate, centroid));

    if (clockwise == true)
    {
        double order = fixPlaneNormalIn.dot((points[0]).cross(points[1]));

        if (order > 0)
        {
            std::reverse(points.begin(), points.end());
        }
    }
}

void Plane::reverse()
{
    normalVector = -1.0 * normalVector;
    bias = -1.0 * bias;
}

void Plane::setCenter(const Point& pPoint)
{
	mPoint = pPoint;
	this->bias = (-1.0) * normalVector.dot(pPoint);
}

void Plane::reverseByNormal(const Point& checkNormal)
{
    if (checkNormal.dot(normalVector) < 0)
    {
        reverse();
    }
}

void Plane::reverseByPoint(const Point& pPoint, bool sameDirection)
{
    if (sameDirection == true)
    {
        if (eval(pPoint) < 0)
        {
            reverse();
        }
    }
    else
    {
        if (eval(pPoint) > 0)
        {
            reverse();
        }
    }
}

void Plane::movePlaneOnNormal(double distance)
{
    Point newPoint = mPoint + distance * normalVector;
    mPoint = newPoint;
    bias = (-1.0)*normalVector.dot(newPoint);
}

double Plane::eval(const double a[3]) const
{
    Point myPoint(a[0], a[1], a[2]);
    return eval(myPoint);
}

void Plane::countPositiveAndNegativePoints(const vtkSmartPointer<vtkPoints>& vtkPointsList, int& positive, int& negative) const
{
	positive = 0;
	negative = 0;

	int tSize = vtkPointsList->GetNumberOfPoints();

	for (int i = 0; i < tSize; i++)
	{
		double pnt[3];
		vtkPointsList->GetPoint(i, pnt);
		if (eval(pnt) >= 0)
		{
			positive++;
		}
		else
		{
			negative++;
		}
	}
}

void Plane::show() const
{
	std::cout << "Normal: " << std::setprecision(15) << normalVector.x << ", " << normalVector.y << ", " << normalVector.z << std::endl;
	std::cout << "Bias: " << std::setprecision(15) << bias << std::endl;
	std::cout << "\n" << std::endl;
}
