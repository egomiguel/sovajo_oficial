#ifndef SLICE_BORDER_PRIVATE_H
#define SLICE_BORDER_PRIVATE_H

#include <opencv2/calib3d/calib3d.hpp>
#include "Segment.hpp"
#include <unsupported/Eigen/Splines>
#include <list>

namespace TKA
{
	namespace SEGMENTATION
	{
		class SliceBorderPrivate
		{
		public:
			SliceBorderPrivate();
			bool DeletePointsInsideRadius(std::list<cv::Point3d>& points, const cv::Point3d& centerPoint, cv::Point3d& nearPoint, double radius);
			//std::vector<cv::Point3d> SortPoints(const std::list<cv::Point3d>& points);
			int GetLessAnglePosition(const std::vector<JoinedSegment>& segments, double& angle);
			std::vector<cv::Point3d> GetContourMainPoints(const std::vector<cv::Point3d>& sortPoints, double angle);
			std::vector<cv::Point3d> InterpolateSpline(const std::vector<cv::Point3d>& pPoints, int amount);
			double GetDistanceBetweenPoints(cv::Point3d a, cv::Point3d b);
			cv::Point3d ArrayToPoint(const double point[3]);
		};
	}
}

#endif