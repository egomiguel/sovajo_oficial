#ifndef THA_HIP_CENTER_H
#define THA_HIP_CENTER_H

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
			static Sphere TestHipCenterBySphere(const cv::Point3d& center, double radius, bool addNoise);
			
		private:
			std::vector<std::vector<cv::Point3d> > ellipse_list_;
			std::vector<cv::Point3d> point_list_;
			bool is_ellipse_;
			static cv::Mat LeastSquareSolve(const cv::Mat& A, const cv::Mat& B);
			static Sphere GetSphereCenter(const std::vector<cv::Point3d>& sphere_points);
			static std::vector<double> computeResiduals(const std::vector<cv::Point3d>& points, const cv::Point3d& center, double radius);
			static double computeMAD(const std::vector<double>& residuals);
			static double median(std::vector<double> values);
			static std::vector<cv::Point3d> removeOutliersMAD(
				const std::vector<cv::Point3d>& points,
				const cv::Point3d& center,
				double radius,
				double zThreshold = 1.0,
				double maxRemovalFraction = 0.05);
			static Sphere RefineSphere(const cv::Point3d& initialCenter, double initialRadius, const std::vector<cv::Point3d>& points, int maxIterations);
			static std::pair<std::vector<cv::Point3d>, double> GenerateHemisphere(const cv::Point3d& center, double radius, int numPoints, bool addNoise);
		};

	}
}

#endif
