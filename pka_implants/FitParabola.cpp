#include "FitParabola.hpp"
#include "Line.hpp"

using namespace PKA::IMPLANTS;

FitParabola::FitParabola(const std::vector<Point>& pPoints, const Plane& plane, bool& result)
{
    Point zAxis(0, 0, 1);
    Point rotAxis = (plane.getNormalVector()).cross(zAxis);
    double angle = Line::getAngleBetweenVectors(zAxis, plane.getNormalVector());

    cv::Mat rotation_1 = Line::getRotateMatrix(rotAxis, -angle);
    cv::Mat rotation_2 = Line::getRotateMatrix(rotAxis, angle);
    cv::Mat rotateVector_1 = rotation_1 * plane.getNormalVectorMat();
    cv::Mat rotateVector_2 = rotation_2 * plane.getNormalVectorMat();

    double distance_1 = Line::getAngleBetweenVectors(Point(rotateVector_1), zAxis);
    double distance_2 = Line::getAngleBetweenVectors(Point(rotateVector_2), zAxis);

    if (distance_1 < distance_2)
    {
        zRotation = rotation_1;
    }
    else
    {
        zRotation = rotation_2;
    }

    zCoord = 0;

    for (int i = 0; i < pPoints.size(); i++)
    {
        cv::Mat pointMat = zRotation * (pPoints[i].ToMatPoint());
        Point myPoint = Point(pointMat);
        points.push_back(myPoint);
        zCoord += myPoint.z;
    }

    if (pPoints.size() > 0)
    {
        zCoord = zCoord / double(pPoints.size());
    }

    result = makeFit();
}

bool FitParabola::makeFit()
{
    double x = 0, y = 0, x2 = 0, x3 = 0, x4 = 0, xy = 0, x2y = 0;

    for (int i = 0; i < points.size(); i++)
    {
        x += points[i].x;
        y += points[i].y;
        x2 += points[i].x * points[i].x;
        x3 += points[i].x * points[i].x * points[i].x;
        x4 += points[i].x * points[i].x * points[i].x * points[i].x;
        xy += points[i].x * points[i].y;
        x2y += points[i].x * points[i].x * points[i].y;
    }

    cv::Mat A(3, 3, CV_64F);
    cv::Mat B(3, 1, CV_64F);
    cv::Mat S;

    A.at<double>(0, 0) = double(points.size());
    A.at<double>(0, 1) = x;
    A.at<double>(0, 2) = x2;

    A.at<double>(1, 0) = x;
    A.at<double>(1, 1) = x2;
    A.at<double>(1, 2) = x3;

    A.at<double>(2, 0) = x2;
    A.at<double>(2, 1) = x3;
    A.at<double>(2, 2) = x4;

    B.at<double>(0, 0) = y;
    B.at<double>(1, 0) = xy;
    B.at<double>(2, 0) = x2y;

    bool result = cv::solve(A, B, S);
    Point pointResult = Point(S);
    a = pointResult.z;
    b = pointResult.y;
    c = pointResult.x;

    if (abs(a) < 0.0000001)
    {
        result = false;
    }

    return result;
}

Point FitParabola::getVertex() const
{
    double Xv = -b / (2.0 * a);
    double Yv = a * Xv * Xv + b * Xv + c;
    Point Pv(Xv, Yv, zCoord);
    cv::Mat resultMat = zRotation.inv() * Pv.ToMatPoint();
    return Point(resultMat);
}

std::vector<Point> FitParabola::getPoints(double rangeInit, double rangeEnd) const
{
    std::vector<Point> result;

    for (double i = rangeInit; i <= rangeEnd; i++)
    {
        double y = a * i * i + b * i + c;
        Point Pv(i, y, zCoord);
        cv::Mat resultMat = zRotation.inv() * Pv.ToMatPoint();
        result.push_back(Point(resultMat));
    }
    result.push_back(getVertex());
    return result;
}

std::vector<Point> FitParabola::SmoothCurve(const std::vector<Point>& pPoints, int kernelSize)
{
    std::vector<Point> result;
    int tSize = pPoints.size();

    if (tSize == 0 || kernelSize == 0)
    {
        return result;
    }

    int kernel;
    if ((tSize / 5) >= kernelSize)
    {
        kernel = kernelSize;
    }
    else
    {
        kernel = tSize / 5;
    }
    //std::cout << "Kernel: " << kernel << std::endl;

    for (int i = 0; i < tSize - kernel - 1; i++)
    {
        Point suma;
        for (int j = i; j < i + kernel; j++)
        {
            suma = suma + pPoints[j];
        }
        suma = suma / double(kernel);
        result.push_back(suma);
    }

    //result.push_back(pPoints[pPoints.size() - 1]);
    return result;
}