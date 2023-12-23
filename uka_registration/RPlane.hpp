#ifndef REGISTRATION_PLANE_H
#define REGISTRATION_PLANE_H

#include "uka_registration_export.h"
#include "RLine.hpp"

namespace UKA
{
	namespace REGISTRATION
	{
		class UKA_REGISTRATION_EXPORT RPlane
		{
			//ax + by + cz + d = 0
		public:
			RPlane();

			void init(const cv::Point3d& point1, const cv::Point3d& point2, const cv::Point3d& point3, bool normalice = false);

			void init(const cv::Point3d& normalVector, const cv::Point3d& pPoint, bool normalice = false);

			double getSquareNorm(const cv::Point3d& pPoint) const;

			void normalizeNormalVector();

			cv::Mat getNormalVectorMat() const;

			cv::Mat getPointMat() const;

			RPlane getPerpendicularPlane(const cv::Point3d& point1, const cv::Point3d& point2) const;

			cv::Point3d getInterceptionLinePoint(const RLine& a) const;

			cv::Point3d getProjectionVector(const cv::Point3d& vector) const;

			cv::Point3d getProjectionPoint(const cv::Point3d& pPoint) const;

			bool isPointBelongToPlane(const cv::Point3d& pPoint) const;

			bool isPointNearToPlane(const cv::Point3d& pPoint, double distance = 0.5) const;

			double getDistanceFromPoint(const cv::Point3d& pPoint) const;

			cv::Point3d getNormalVector() const;

			double getBias() const;

			cv::Point3d getPoint() const;

			static RPlane getBestPlane(const std::vector<cv::Point3d>& pPoints);

			double eval(const cv::Point3d& a) const;

			double eval(const double a[3]) const;

			void show() const;

			void reverse();

			void reverseNormal(const cv::Point3d& a, bool sameDirection);

			bool getIsInit() const;

		private:
			cv::Point3d normalVector;
			cv::Point3d mPoint;
			double bias;
			bool isInit;
			void fixNormalVector(cv::Point3d& newNormalVector);
			void setPoint(cv::Point3d& pPoint);
		};
	}
}

#endif