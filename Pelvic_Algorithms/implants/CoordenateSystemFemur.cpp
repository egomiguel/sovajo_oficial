#include "CoordenateSystemFemur.hpp"

CoordenateSystemFemur::CoordenateSystemFemur(const Point& hipCenter, const Point& latEpicondyle, const Point& medEpicondyle, const Point& kneeCenter)
{
    Point vectorZ = hipCenter - kneeCenter;
    Point vectorX = latEpicondyle - medEpicondyle;
    vectorZ.normalice();
    vectorX.normalice();
    Point vectorY = vectorX.cross(vectorZ);
    X.init(vectorX, kneeCenter);
    Y.init(vectorY, kneeCenter);
    Z.init(vectorZ, kneeCenter);
}

Point CoordenateSystemFemur::getPointCoordenate(const Point& pPoint) const
{
    double a = X.eval(pPoint);
    double b = Y.eval(pPoint);
    double c = Z.eval(pPoint);
    return Point(a, b, c);
}

Plane CoordenateSystemFemur::getTransformPlane(const Plane& myPlane, const cv::Mat& rotation, const cv::Mat& translation) const
{
    cv::Mat normalMat = myPlane.getNormalVectorMat() * rotation;
    cv::Mat pointMat = myPlane.getPointMat() * rotation + translation;
    Plane result;
    result.init(Point(normalMat), Point(pointMat));
    return result;
}

Point CoordenateSystemFemur::getPointCoordenate(const Point& pPoint, const cv::Mat& rotation, const cv::Mat& translation) const
{
    Plane XT = getTransformPlane(X, rotation, translation);
    Plane YT = getTransformPlane(Y, rotation, translation);
    Plane ZT = getTransformPlane(Z, rotation, translation);

    cv::Mat pointMat = rotation * (pPoint.ToMatPoint()) + translation;
    Point newPoint = Point(pointMat);

    double a = XT.eval(newPoint);
    double b = YT.eval(newPoint);
    double c = ZT.eval(newPoint);

    return Point(a, b, c);
}

Point CoordenateSystemFemur::getPointCoordenate(const Point& pPoint, const cv::Mat& planeTransform, bool transformPoint) const
{
    double* matrix = new double[9];

    matrix[0] = planeTransform.at<double>(0, 0);
    matrix[1] = planeTransform.at<double>(0, 1);
    matrix[2] = planeTransform.at<double>(0, 2);

    matrix[3] = planeTransform.at<double>(1, 0);
    matrix[4] = planeTransform.at<double>(1, 1);
    matrix[5] = planeTransform.at<double>(1, 2);

    matrix[6] = planeTransform.at<double>(2, 0);
    matrix[7] = planeTransform.at<double>(2, 1);
    matrix[8] = planeTransform.at<double>(2, 2);

    cv::Mat rotation(3, 3, CV_64FC1, matrix);

    cv::Mat translation(3, 1, CV_64FC1);
    translation.at <double>(0, 0) = planeTransform.at<double>(0, 3);
    translation.at <double>(1, 0) = planeTransform.at<double>(1, 3);
    translation.at <double>(2, 0) = planeTransform.at<double>(2, 3);

    Plane XT = getTransformPlane(X, rotation, translation);
    Plane YT = getTransformPlane(Y, rotation, translation);
    Plane ZT = getTransformPlane(Z, rotation, translation);

    double a, b, c;

    if (transformPoint == true)
    {
        cv::Mat pointMat = rotation * (pPoint.ToMatPoint()) + translation;
        Point newPoint = Point(pointMat);

        a = XT.eval(newPoint);
        b = YT.eval(newPoint);
        c = ZT.eval(newPoint);
    }
    else
    {
        a = XT.eval(pPoint);
        b = YT.eval(pPoint);
        c = ZT.eval(pPoint);
    }

    return Point(a, b, c);
}

void CoordenateSystemFemur::getPointCoordenate(const std::vector<Point>& pointsIn, std::vector<PointTypeITK>& pointsOut, const cv::Mat& rotation, const cv::Mat& translation) const
{
    Plane XT = getTransformPlane(X, rotation, translation);
    Plane YT = getTransformPlane(Y, rotation, translation);
    Plane ZT = getTransformPlane(Z, rotation, translation);

    double a, b, c;
    for (int i = 0; i < pointsIn.size(); i++)
    {

        cv::Mat pointMat = rotation * (pointsIn[i].ToMatPoint()) + translation;
        Point newPoint = Point(pointMat);

        a = XT.eval(newPoint);
        b = YT.eval(newPoint);
        c = ZT.eval(newPoint);

        PointTypeITK myPoint;
        myPoint[0] = a;
        myPoint[1] = b;
        myPoint[2] = c;

        pointsOut.push_back(myPoint);
    }
}

void CoordenateSystemFemur::getPointCoordenate(const std::vector<Point>& pointsIn, std::vector<PointTypeITK>& pointsOut, const cv::Mat& planeTransform) const
{
    double* matrix = new double[9];

    matrix[0] = planeTransform.at<double>(0, 0);
    matrix[1] = planeTransform.at<double>(0, 1);
    matrix[2] = planeTransform.at<double>(0, 2);

    matrix[3] = planeTransform.at<double>(1, 0);
    matrix[4] = planeTransform.at<double>(1, 1);
    matrix[5] = planeTransform.at<double>(1, 2);

    matrix[6] = planeTransform.at<double>(2, 0);
    matrix[7] = planeTransform.at<double>(2, 1);
    matrix[8] = planeTransform.at<double>(2, 2);

    cv::Mat rotation(3, 3, CV_64FC1, matrix);

    cv::Mat translation(3, 1, CV_64FC1);
    translation.at <double>(0, 0) = planeTransform.at<double>(0, 3);
    translation.at <double>(1, 0) = planeTransform.at<double>(1, 3);
    translation.at <double>(2, 0) = planeTransform.at<double>(2, 3);

    Plane XT = getTransformPlane(X, rotation, translation);
    Plane YT = getTransformPlane(Y, rotation, translation);
    Plane ZT = getTransformPlane(Z, rotation, translation);

    double a, b, c;
    for (int i = 0; i < pointsIn.size(); i++)
    {

        cv::Mat pointMat = rotation * (pointsIn[i].ToMatPoint()) + translation;
        Point newPoint = Point(pointMat);

        a = XT.eval(newPoint);
        b = YT.eval(newPoint);
        c = ZT.eval(newPoint);

        PointTypeITK myPoint;
        myPoint[0] = a;
        myPoint[1] = b;
        myPoint[2] = c;

        pointsOut.push_back(myPoint);
    }
}