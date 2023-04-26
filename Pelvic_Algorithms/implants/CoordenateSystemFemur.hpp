#include "Point.hpp"
#include "Plane.hpp"

class CoordenateSystemFemur
{
private:
    Plane X, Y, Z;
    Plane getTransformPlane(const Plane& myPlane, const cv::Mat& rotation, const cv::Mat& translation) const;
public:
    CoordenateSystemFemur(const Point& hipCenter, const Point& latEpicondyle, const Point& medEpicondyle, const Point& kneeCenter);
    Point getPointCoordenate(const Point& pPoint) const;
    Point getPointCoordenate(const Point& pPoint, const cv::Mat& planeTransform, bool transformPoint = true) const;
    Point getPointCoordenate(const Point& pPoint, const cv::Mat& rotation, const cv::Mat& translation) const;
    void getPointCoordenate(const std::vector<Point>& pointsIn, std::vector<PointTypeITK>& pointsOut, const cv::Mat& planeTransform) const;
    void getPointCoordenate(const std::vector<Point>& pointsIn, std::vector<PointTypeITK>& pointsOut, const cv::Mat& rotation, const cv::Mat& translation) const;
};