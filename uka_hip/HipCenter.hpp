#ifndef HIP_CENTER_H
#define HIP_CENTER_H

#include "DataFit.hpp"
#include "uka_hip_export.h"

namespace UKA
{
	namespace HIP
	{

		struct HipPoint
		{
			double x, y, z;
		};

		class UKA_HIP_EXPORT HipCenter
		{
		public:
			HipCenter(const std::vector<std::vector<HipPoint> >& ellipse_list);
			HipCenter(const std::vector<HipPoint>& point_list);
			void GetHipCenterByEllipses(double point_out[3]) const;
			void GetHipCenterBySphere(double point_out[3]) const;
			void GetHipCenterByThreeSphere(double point_out[3]) const;
			static void GetHipCenterFromFile(std::string path, double point_out[3]);

		private:
			std::vector<std::vector<cv::Point3d> > ellipse_list_;
			std::vector<cv::Point3d> point_list_;
			bool is_ellipse_;
			cv::Point3d LeastSquareSolve(const cv::Mat& A, const cv::Mat& B) const;
			cv::Point3d GetSphereCenter(const std::vector<cv::Point3d>& sphere_points) const;
		};

	}
}

#endif
