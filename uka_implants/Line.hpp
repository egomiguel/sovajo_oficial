#ifndef LINE_UKA_H
#define LINE_UKA_H

#include "Point.hpp"
#include "uka_implants_export.h"

namespace UKA
{
	namespace IMPLANTS
	{
		class UKA_IMPLANTS_EXPORT Line
		{
		public:
			Line(const Point& directVector, const Point& pPoint);

			Line(const Line& pLine);

			Point getPoint() const;

			Point getPointAtDistance(double distance) const;

			Point getDirectVector() const;

			double getSquareNorm(const Point& pPoint) const;

			double getDistanceFromPoint(const Point& pPoint) const;

			Line getParalellLine(const Point& pPoint) const;

			Line getPerpendicularLine(const Point& pPoint) const;

			bool isPointBelongToLine(const Point& pPoint) const;

			Point getPointAtDistance(const Point& pPoint, const Point& nearReferencePoint, float distance, bool closest = true) const;

			static Point getFixDirectVector(const Point& pLineDirectVector, const Point& pLineFixPoint, const Point& pReferencePoint, bool closest = true);

			static cv::Mat getRotateMatrix(const Point& axis, double angle);

			static double getAngleBetweenVectors(const Point& a, const Point& b);

			static double getDistanceBetweenPoints(const Point& a, const Point& b, bool square = false);

			static Point getProjectPoint(const Point& linePoint1, const Point& linePoint2, const Point& externalPoint);

			Point getProjectPoint(const Point& pPoint) const;

			void setPoint(const Point& newPoint);

			void setDirectVector(const Point& newVector);

			void normaliceDirectVector();

			bool getInterceptionSphere(const Point& center, double radius, std::pair<Point, Point>& result) const;

			Point getNearestPoint(const std::vector<Point>& points);

			Point evalParameter(double t) const;

			static Line makeLineWithPoints(const Point& P1, const Point& P2);

		private:
			Point mPoint;
			Point directVector;
		};
	}
}

#endif
