#ifndef REGISTRATION_LINE_H
#define REGISTRATION_LINE_H

#include <opencv2/calib3d/calib3d.hpp>
#include "uka_registration_export.h"

namespace UKA
{
	namespace REGISTRATION
	{
		class UKA_REGISTRATION_EXPORT RLine
		{
		public:
			RLine(const cv::Point3d& directVector, const cv::Point3d& pPoint);

			cv::Point3d getPoint() const;

			cv::Point3d getDirectVector() const;

			double getSquareNorm(const cv::Point3d& pPoint) const;

			double getDistanceFromPoint(const cv::Point3d& pPoint) const;

			double getDistanceFromPoint(const double pPoint[3]) const;

			RLine getParalellLine(const cv::Point3d& pPoint) const;

			RLine getPerpendicularLine(const cv::Point3d& pPoint) const;

			cv::Point3d getPointAtDistance(const cv::Point3d& pPoint, const cv::Point3d& nearReferencePoint, float distance, bool closest = true) const;

			cv::Point3d getProjectPoint(const cv::Point3d& pPoint) const;

			void normaliceDirectVector();

			static double getDistanceBetweenPoints(const cv::Point3d& a, const cv::Point3d& b, bool square = false);

			bool isPointBelongToLine(const cv::Point3d& pPoint) const;

			static cv::Mat getRotateMatrix(const cv::Point3d& axis, double angle);

			static double getAngleBetweenVectors(const cv::Point3d& a, const cv::Point3d& b);

			static RLine makeLineWithPoints(const cv::Point3d& P1, const cv::Point3d& P2);

			bool getInterceptionSphere(const cv::Point3d& center, double radius, std::pair<cv::Point3d, cv::Point3d>& result) const;

			cv::Point3d evalParameter(double t) const;

		private:
			cv::Point3d mPoint;
			cv::Point3d directVector;
		};
	}
}

#endif
