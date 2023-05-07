#include <iostream>
#include <math.h>
#include <iomanip>
#include "SPlane.hpp"
#include "SegmentationException.hpp"

using namespace TKA::SEGMENTATION;

const double EPSILON = std::numeric_limits<double>::epsilon();

SPlane::SPlane()
{
    normalVector = cv::Point3d(0.0, 0.0, 0.0);
    bias = 0.0;
    isInit = false;
}

bool SPlane::GetIsInit() const
{
    return isInit;
}

void SPlane::init(const double point1[3], const double point2[3], const double point3[3])
{
    cv::Point3d P1, P2, P3;
    P1 = cv::Point3d(point1[0], point1[1], point1[2]);
    P2 = cv::Point3d(point2[0], point2[1], point2[2]);
    P3 = cv::Point3d(point3[0], point3[1], point3[2]);
    init(P1, P2, P3);
}

void SPlane::init(const double normalVector[3], const double pPoint[3])
{
    cv::Point3d P1, P2;
    P1 = cv::Point3d(normalVector[0], normalVector[1], normalVector[2]);
    P2 = cv::Point3d(pPoint[0], pPoint[1], pPoint[2]);
    init(P1, P2);
}

void SPlane::init(const cv::Point3d& point1, const cv::Point3d& point2, const cv::Point3d& point3)
{
    cv::Point3d directVector1, directVector2;
    directVector1 = point2 - point1;
    directVector2 = point3 - point1;

    cv::Point3d crossProduct = directVector1.cross(directVector2);

    if (abs(crossProduct.dot(crossProduct)) < EPSILON)
    {
        throw SegmentationExceptionCode::POINTS_TO_DEFINE_PLANE_CAN_NOT_BE_COLLINEAR_ON_SEGMENTATION;
    }

    if (isInit == true)
    {
        throw SegmentationExceptionCode::ALREADY_INITIALIZED_PLANE_ON_SEGMENTATION;
    }

    mPoint = point3;
    normalVector = directVector1.cross(directVector2);
    bias = (-1.0)*normalVector.dot(point3);
    normalizeNormalVector();
    isInit = true;
}

void SPlane::init(const cv::Point3d& normalVector, const cv::Point3d& pPoint)
{
    if (isInit == true)
    {
        throw SegmentationExceptionCode::ALREADY_INITIALIZED_PLANE_ON_SEGMENTATION;
    }

    if (normalVector.dot(normalVector) == 0)
    {
        throw SegmentationExceptionCode::PLANE_NORMAL_VECTOR_CAN_NOT_BE_ZERO_ON_SEGMENTATION;
    }

    mPoint = pPoint;
    this->normalVector = normalVector;
    this->bias = (-1) * normalVector.dot(pPoint);
    normalizeNormalVector();
    isInit = true;
}

void SPlane::normalizeNormalVector()
{
    double squareNorm = sqrt(normalVector.dot(normalVector));
    normalVector = (normalVector / squareNorm);
    bias = (bias / squareNorm);
}


SPlane SPlane::getPerpendicularPlane(const cv::Point3d& point1, const cv::Point3d& point2) const
{
    cv::Point3d tDirectVector = point2 - point1;
    cv::Point3d tNormalVector = normalVector.cross(tDirectVector);

    if (tNormalVector.dot(tNormalVector) < EPSILON)
    {
        throw SegmentationExceptionCode::CAN_NOT_DETERMINE_PERPENDICULAR_PLANE_TO_PLANE_ON_SEGMENTATION;
    }

    SPlane perpendicularPlane;
    perpendicularPlane.init(tNormalVector, point1);
    return perpendicularPlane;
}

SPlane SPlane::getPerpendicularPlane(const double point1[3], const double point2[3]) const
{
    cv::Point3d newP1 = cv::Point3d(point1[0], point1[1], point1[2]);
    cv::Point3d newP2 = cv::Point3d(point2[0], point2[1], point2[2]);

    return getPerpendicularPlane(newP1, newP2);
}

cv::Point3d SPlane::getProjectionPoint(const cv::Point3d& pPoint) const
{
    cv::Point3d diff = pPoint - mPoint;
    cv::Point3d projectOnNormal = ((diff.dot(normalVector)) / (normalVector.dot(normalVector))) * normalVector;
    cv::Point3d projection = pPoint - projectOnNormal;
    return projection;
}

cv::Point3d SPlane::getProjectionPoint(const double p[3]) const
{
    cv::Point3d point = cv::Point3d(p[0], p[1], p[2]);
    cv::Point3d proj = getProjectionPoint(point);
    return proj;
}

cv::Point3d SPlane::getProjectionVector(const cv::Point3d& vector) const
{
    cv::Point3d projectOnNormal = ((vector.dot(normalVector)) / (normalVector.dot(normalVector))) * normalVector;
    cv::Point3d projection = vector - projectOnNormal;
    return projection;
}

double SPlane::getDistanceFromPoint(const cv::Point3d& pPoint) const
{
    double squareNorm = normalVector.dot(normalVector);
    double distance = abs(normalVector.dot(pPoint) + bias) / sqrt(squareNorm);
    return distance;
}

cv::Point3d SPlane::getNormalVector() const
{
    return normalVector;
}

cv::Point3d SPlane::getPoint() const
{
    return mPoint;
}

double SPlane::eval(const cv::Point3d& p) const
{
    return normalVector.dot(p) + bias;
}

double SPlane::eval(const double p[3]) const
{
    cv::Point3d point = cv::Point3d(p[0], p[1], p[2]);

    return normalVector.dot(point) + bias;
}

void SPlane::reverse()
{
    normalVector = -1.0 * normalVector;
    bias = -1.0 * bias;
}

bool SPlane::isParallelToPlane(const SPlane& plane) const
{
    cv::Point3d externalNormal = plane.getNormalVector();

    cv::Point3d crossProduct = externalNormal.cross(normalVector);

    if (abs(crossProduct.dot(crossProduct)) < EPSILON)
    {
        return true;
    }
    return false;
}

void SPlane::move(const double p[3])
{
    cv::Point3d point = cv::Point3d(p[0], p[1], p[2]);
    mPoint = point;
    bias = (-1.0)*normalVector.dot(point);
}

void SPlane::move(const cv::Point3d& p)
{
    mPoint = p;
    bias = (-1.0)*normalVector.dot(p);
}

