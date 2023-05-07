#ifndef POINT_H
#define POINT_H

#include <opencv2/calib3d/calib3d.hpp>
#include "tka_implants_export.h"
#include "Types.hpp"
namespace TKA
{
	namespace IMPLANTS
	{
		class TKA_IMPLANTS_EXPORT Point : public cv::Point3d
		{
		public:
			Point();
			Point(double x, double y, double z);
			Point(const cv::Point3d& P);
			Point(const Point& P);
			Point(const cv::Mat& M);
			Point& operator = (const cv::Point3d& P);
			Point& operator = (const Point& P);
			cv::Point3d ToCVPoint() const;
			cv::Mat ToMatPoint() const;
			PointTypeITK ToITKPoint() const;
			void normalice();
		};
	}
}
#endif