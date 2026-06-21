#ifndef HIP_CENTER_H
#define HIP_CENTER_H

#include "DataFit.hpp"
#include "tha_hip_export.h"

namespace THA
{
	namespace HIP
	{

		struct HipPoint
		{
			double x, y, z;
		};

		class THA_HIP_EXPORT HipCenter
		{
		public:
			struct Sphere
			{
				cv::Point3d center;
				double radius;
				double error;
			};

			HipCenter(const std::vector<std::vector<HipPoint> >& ellipse_list);
			HipCenter(const std::vector<HipPoint>& point_list);
			void GetHipCenterByEllipses(double point_out[3]) const;
			Sphere GetHipCenterBySphere() const;
			Sphere GetHipCenterFromFile(std::string path);
			static std::vector<double> computeResiduals(const std::vector<cv::Point3d>& points, const cv::Point3d& center, double radius);
			static double computeMAD(const std::vector<double>& residuals);
			static double median(std::vector<double> values);
			static std::vector<cv::Point3d> removeOutliersMAD(
				const std::vector<cv::Point3d>& points,
				const cv::Point3d& center,
				double radius,
				double zThreshold = 3.0,
				double maxRemovalFraction = 0.05);
			static Sphere RefineSphere(const cv::Point3d& initialCenter, double initialRadius, const std::vector<cv::Point3d>& points, int maxIterations);

		private:
			std::vector<std::vector<cv::Point3d> > ellipse_list_;
			std::vector<cv::Point3d> point_list_;
			bool is_ellipse_;
			cv::Mat LeastSquareSolve(const cv::Mat& A, const cv::Mat& B) const;
			Sphere GetSphereCenter(const std::vector<cv::Point3d>& sphere_points) const;
			
		};

	}
}

#endif
