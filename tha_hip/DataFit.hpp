#ifndef DATA_FIT_H
#define DATA_FIT_H

#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>

namespace THA
{
	namespace HIP
	{

		struct Normal
		{
			cv::Point3d equation1;
			cv::Point3d equation2;
			double bias1, bias2;
		};

		class DataFit
		{
			friend class HipCenter;
		private:
			std::vector<cv::Point3d> ellipse_points_;
			cv::Point3d centroid_;
			std::vector<double> plane_;
			Normal normal_;
			cv::Point3d GetCentroid(const std::vector<cv::Point3d>& points);
			std::vector<double> FitPlane(std::vector<cv::Point3d>& points, const cv::Point3d& centroid);
			Normal GetNormalEquation();
		public:
			DataFit(const std::vector<cv::Point3d>& ellipse_points);
		};
	}
}

#endif