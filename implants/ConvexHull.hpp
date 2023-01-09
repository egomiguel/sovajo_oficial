#ifndef REGION2D_IMPLANT_H
#define REGION2D_IMPLANT_H

#include <opencv2/calib3d/calib3d.hpp>
#include "Segment.hpp"
#include "Point.hpp"


class ConvexHull
{
private:
    std::vector<cv::Point2f> points;
    cv::Mat rotationOnZInv;
    float axisZ;
    //std::vector<Point> GetPoints3D(const std::vector<cv::Point2f>& pPoints);
public:
    ConvexHull(const std::vector<Point>& pPoints, const cv::Mat& pRotationOnZ);

    std::vector<Point> GetConvexHull(int vertices = -1);

    static std::vector<Point> ReduceConvexHull(const std::vector<Point>& pPoints, int vertices = -1);

    static int GetLessAnglePosition(const std::vector<JoinedSegment>& segments);

    static std::vector<Point> interpolateSpline(const std::vector<Point>& pPoints, int amount = 200);
};


#endif