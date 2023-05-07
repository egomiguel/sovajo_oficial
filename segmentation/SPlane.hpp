#ifndef SEGMENTATION_PLANE_H
#define SEGMENTATION_PLANE_H

#include <opencv2/calib3d/calib3d.hpp>
#include "tka_segmentation_export.h"

namespace TKA
{
	namespace SEGMENTATION
	{
		class TKA_SEGMENTATION_EXPORT SPlane
		{
			//ax + by + cz + d = 0
		public:
			SPlane();

			void init(const cv::Point3d& point1, const cv::Point3d& point2, const cv::Point3d& point3);

			void init(const cv::Point3d& normalVector, const cv::Point3d& pPoint);

			void init(const double point1[3], const double point2[3], const double point3[3]);

			void init(const double normalVector[3], const double pPoint[3]);

			cv::Point3d getPoint() const;

			SPlane getPerpendicularPlane(const cv::Point3d& point1, const cv::Point3d& point2) const;

			SPlane getPerpendicularPlane(const double point1[3], const double point2[3]) const;

			bool isParallelToPlane(const SPlane& plane) const;

			cv::Point3d getProjectionPoint(const cv::Point3d& pPoint) const;

			cv::Point3d getProjectionPoint(const double p[3]) const;

			cv::Point3d getProjectionVector(const cv::Point3d& vector) const;

			cv::Point3d getNormalVector() const;

			double getDistanceFromPoint(const cv::Point3d& pPoint) const;

			double eval(const cv::Point3d& p) const;

			double eval(const double p[3]) const;

			void reverse();

			void move(const cv::Point3d& p);

			void move(const double p[3]);

			bool GetIsInit() const;

		private:
			cv::Point3d normalVector;
			cv::Point3d mPoint;
			double bias;
			bool isInit;
			void normalizeNormalVector();
		};
	}
}

#endif