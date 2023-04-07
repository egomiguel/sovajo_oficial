#ifndef REGION2D_IMPLANT_H
#define REGION2D_IMPLANT_H

#include <opencv2/calib3d/calib3d.hpp>
#include "Segment.hpp"
#include "Point.hpp"


class ConvexHull
{
private:
    std::vector<cv::Point2f> points;
	std::vector<cv::Point2f> mConvexHull2D;
	std::vector<Point> mConvexHull;
    cv::Mat rotationOnZ, rotationOnZInv;
	cv::Point2f mCenterPoint;
    float axisZ;
	std::vector<cv::Point2f> increaseConvexHull(float increaseDist);
    //std::vector<Point> GetPoints3D(const std::vector<cv::Point2f>& pPoints);
public:
    ConvexHull(const std::vector<Point>& pPoints, const cv::Mat& pRotationOnZ);

	bool isPointWithinConvexHull(const Point& pPoint, float pErrorMargin = 0);

	bool areSomePointWithinConvexHull(const std::vector<Point>& pPoints, float pErrorMargin = 0);

	std::vector<Point> getChangeDirectionPoints();

    std::vector<Point> GetConvexHull(int vertices = -1);

    static std::vector<Point> ReduceConvexHull(const std::vector<Point>& pPoints, int vertices = -1);

    static int GetLessAnglePosition(const std::vector<JoinedSegment>& segments);

    static std::vector<Point> interpolateSpline(const std::vector<Point>& pPoints, int amount = 200);
};


#endif