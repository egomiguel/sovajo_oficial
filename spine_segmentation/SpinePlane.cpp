#include <iostream>
#include <math.h>
#include <iomanip>
#include "SpinePlane.hpp"
#include "SegmentationException.hpp"

using namespace SPINE::SEGMENTATION;

const double EPSILON = std::numeric_limits<double>::epsilon();

SpinePlane::SpinePlane()
{
    normalVector = cv::Point3d(0.0, 0.0, 0.0);
    bias = 0.0;
    isInit = false;
}

void SpinePlane::fixNormalVector(cv::Point3d& newNormalVector)
{
    normalVector = newNormalVector;
    bias = (-1.0)*newNormalVector.dot(mPoint);
    normalizeNormalVector();
}

void SpinePlane::setPoint(cv::Point3d& pPoint)
{
    mPoint = pPoint;
}

void SpinePlane::init(const cv::Point3d& point1, const cv::Point3d& point2, const cv::Point3d& point3, bool normalice)
{
    SpineLine temp((point1 - point2), point2);
    if (temp.isPointBelongToLine(point3))
    {
        throw SegmentationExceptionCode::POINTS_TO_DEFINE_PLANE_CAN_NOT_BE_COLLINEAR;
    }
    cv::Point3d directVector1, directVector2;
    if (isInit == true)
    {
        throw SegmentationExceptionCode::ALREADY_INITIALIZED_PLANE;
    }

    directVector1 = point2 - point1;
    directVector2 = point3 - point1;

    mPoint = point3;
    normalVector = directVector1.cross(directVector2);
    bias = (-1.0)*normalVector.dot(point3);

    normalizeNormalVector();
    isInit = true;
}

void SpinePlane::init(const cv::Point3d& normalVector, const cv::Point3d& pPoint, bool normalice)
{
    if (isInit == true)
    {
        throw SegmentationExceptionCode::ALREADY_INITIALIZED_PLANE;
    }

    mPoint = pPoint;
    this->normalVector = normalVector;
    this->bias = (-1) * normalVector.dot(pPoint);

    normalizeNormalVector();

    isInit = true;
}

double SpinePlane::getSquareNorm(const cv::Point3d& pPoint) const
{
    return pPoint.x * pPoint.x + pPoint.y * pPoint.y + pPoint.z * pPoint.z;
}

void SpinePlane::normalizeNormalVector()
{
    double squareNorm = sqrt(normalVector.dot(normalVector));
    normalVector = (normalVector / squareNorm);
    bias = (bias / squareNorm);
}

cv::Mat SpinePlane::getNormalVectorMat() const
{
    cv::Mat mat(3, 1, CV_64F);
    mat.at <double>(0, 0) = normalVector.x;
    mat.at <double>(1, 0) = normalVector.y;
    mat.at <double>(2, 0) = normalVector.z;
    return mat;
}

cv::Mat SpinePlane::getPointMat() const
{
    cv::Mat mat(3, 1, CV_64FC1);
    mat.at <double>(0, 0) = mPoint.x;
    mat.at <double>(1, 0) = mPoint.y;
    mat.at <double>(2, 0) = mPoint.z;
    return mat;
}

SpinePlane SpinePlane::getPerpendicularPlane(const cv::Point3d& point1, const cv::Point3d& point2) const
{
    cv::Point3d tDirectVector = point2 - point1;
    cv::Point3d tNormalVector = normalVector.cross(tDirectVector);

    if (tNormalVector.dot(tNormalVector) < EPSILON)
    {
        throw SegmentationExceptionCode::CAN_NOT_DETERMINE_PERPENDICULAR_PLANE_TO_PLANE;
    }

    SpinePlane perpendicularPlane;
    perpendicularPlane.init(tNormalVector, point1);
    return perpendicularPlane;
}

cv::Point3d SpinePlane::getInterceptionLinePoint(const SpineLine& a) const
{
    double magnitude = normalVector.dot(a.getDirectVector());
    if (magnitude == 0.0)
    {
        throw SegmentationExceptionCode::CAN_NOT_DETERMINE_INTERCEPTION_LINE_TO_PLANE;
    }
    double parameter = -(bias + (normalVector.dot(a.getPoint()))) / magnitude;

    cv::Point3d result = a.getPoint() + (a.getDirectVector() * parameter);
    return result;
}

cv::Point3d SpinePlane::getProjectionVector(const cv::Point3d& vector) const
{
    cv::Point3d projectOnNormal = ((vector.dot(normalVector)) / (normalVector.dot(normalVector))) * normalVector;
    cv::Point3d projection = vector - projectOnNormal;
    return projection;
}

cv::Point3d SpinePlane::getProjectionPoint(const cv::Point3d& pPoint) const
{
    cv::Point3d diff = pPoint - mPoint;
    cv::Point3d projectOnNormal = ((diff.dot(normalVector)) / (normalVector.dot(normalVector))) * normalVector;
    cv::Point3d projection = pPoint - projectOnNormal;
    return projection;
}

bool SpinePlane::isPointBelongToPlane(const cv::Point3d& pPoint) const
{
    double evaluation = normalVector.dot(pPoint) + bias;
    if (fabs(evaluation) < EPSILON)
    {
        return true;
    }
    return false;
}

bool SpinePlane::isPointNearToPlane(const cv::Point3d& pPoint, double distance) const
{
    double evaluation = (normalVector.dot(pPoint) + bias) / sqrt(getSquareNorm(normalVector));
    if (fabs(evaluation) <= distance)
    {
        return true;
    }
    return false;
}

double SpinePlane::getDistanceFromPoint(const cv::Point3d& pPoint) const
{
    double distance = abs(normalVector.dot(pPoint) + bias) / sqrt(getSquareNorm(normalVector));
    return distance;
}

cv::Point3d SpinePlane::getNormalVector() const
{
    return normalVector;
}

double SpinePlane::getBias() const
{
    return bias;
}

cv::Point3d SpinePlane::getPoint() const
{
    return mPoint;
}

double SpinePlane::eval(const cv::Point3d& a) const
{
    return normalVector.dot(a) + bias;
}

double SpinePlane::eval(const double a[3]) const
{
    cv::Point3d tPoint = cv::Point3d(a[0], a[1], a[2]);
    return eval(tPoint);
}

void SpinePlane::show() const
{
    std::cout << "Normal: " << std::setprecision(15) << normalVector.x << ", " << normalVector.y << ", " << normalVector.z << std::endl;
    std::cout << "Bias: " << std::setprecision(15) << bias << std::endl;
    std::cout << "\n" << std::endl;
}

void SpinePlane::reverse()
{
    normalVector = -1.0 * normalVector;
    bias = -1.0 * bias;
}

void SpinePlane::reverseNormal(const cv::Point3d& a, bool sameDirection)
{
    cv::Point3d proj = getProjectionPoint(a);
    cv::Point3d vector;
    if (sameDirection == true)
    {
        vector = a - proj;
    }
    else
    {
        vector = proj - a;
    }

    fixNormalVector(vector);
}

/*
SpinePlane SpinePlane::getBestPlane(const vtkSmartPointer<vtkPoints> pPoints)
{
    SpinePlane bestPlane;
    int tSize = pPoints->GetNumberOfPoints();
    cv::Point3d averagePoint(0, 0, 0);
    std::vector<cv::Point3d> myPoints;

    for (int i = 0; i < tSize; i++)
    {
        double pnt[3];
        pPoints->GetPoint(i, pnt);
        cv::Point3d temp = { pnt[0], pnt[1], pnt[2] };

        averagePoint = averagePoint + temp;
        myPoints.push_back(temp);
    }

    if (tSize > 0)
    {
        averagePoint = averagePoint / double(tSize);
    }
    else
    {
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

    bestPlane.init(cv::Point3d(a, b, c), averagePoint);
    return bestPlane;
}
*/

SpinePlane SpinePlane::getBestPlane(const std::vector<cv::Point3d>& pPoints)
{
    SpinePlane bestPlane;
    int tSize = pPoints.size();
    cv::Point3d averagePoint(0, 0, 0);
    std::vector<cv::Point3d> myPoints;

    for (int i = 0; i < tSize; i++)
    {
        averagePoint = averagePoint + pPoints[i];
        myPoints.push_back(pPoints[i]);
    }

    if (tSize > 0)
    {
        averagePoint = averagePoint / double(tSize);
    }
    else
    {
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

    bestPlane.init(cv::Point3d(a, b, c), averagePoint);
    return bestPlane;
}

bool SpinePlane::getIsInit() const
{
    return isInit;
}