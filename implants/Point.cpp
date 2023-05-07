
#include "Point.hpp"

using namespace TKA::IMPLANTS;

Point::Point(double x, double y, double z)
	:cv::Point3d(x, y, z)
{
}

Point::Point()
	: cv::Point3d(0, 0, 0)
{
}

Point::Point(const cv::Point3d& P)
	:cv::Point3d(P)
{
}

Point::Point(const Point& P)
	: cv::Point3d(P.x, P.y, P.z)
{
}

Point::Point(const cv::Mat& M)
{
	*this = cv::Point3d(M);
}

Point& Point::operator = (const cv::Point3d& P)
{
	this->x = P.x;
	this->y = P.y;
	this->z = P.z;
	return *this;
}

Point& Point::operator = (const Point& P)
{
	this->x = P.x;
	this->y = P.y;
	this->z = P.z;
	return *this;
}

cv::Point3d Point::ToCVPoint() const
{
	return cv::Point3d(this->x, this->y, this->z);
}

cv::Mat Point::ToMatPoint() const
{
    cv::Mat pointMat(3, 1, CV_64F);
    pointMat.at <double>(0, 0) = this->x;
    pointMat.at <double>(1, 0) = this->y;
    pointMat.at <double>(2, 0) = this->z;
    return pointMat;
}

PointTypeITK Point::ToITKPoint() const
{
    PointTypeITK pointITK;
    pointITK[0] = this->x;
    pointITK[1] = this->y;
    pointITK[2] = this->z;
    return pointITK;
}

void Point::normalice()
{
    double norm = sqrt((*this).dot(*this));
    this->x = this->x / norm;
    this->y = this->y / norm;
    this->z = this->z / norm;
}