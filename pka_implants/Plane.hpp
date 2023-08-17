#ifndef PLANE_H
#define PLANE_H

#include "Line.hpp"
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include "pka_implants_export.h"

namespace PKA
{
	namespace IMPLANTS
	{
		class PKA_IMPLANTS_EXPORT Plane
		{
			//ax + by + cz + d = 0
		public:
			Plane();

			void init(const Point& point1, const Point& point2, const Point& point3);

			void init(const Point& normalVector, const Point& pPoint);

			double getSquareNorm(const Point& pPoint) const;

			void normalizeNormalVector();

			void transformPlane(const cv::Mat& rotation, const cv::Mat& translation);

			cv::Mat getNormalVectorMat() const;

			cv::Mat getPointMat() const;

			Plane getPerpendicularPlane(const Point& point1, const Point& point2) const;

			Point getInterceptionLinePoint(const Line& a) const;

			Point getProjectionVector(const Point& vector) const;

			Point getProjectionPoint(const Point& pPoint) const;

			bool isPointBelongToPlane(const Point& pPoint) const;

			bool isPointNearToPlane(const Point& pPoint, double distance = 0.5) const;

			double getDistanceFromPoint(const Point& pPoint) const;

			void countPositiveAndNegativePoints(const vtkSmartPointer<vtkPoints>& vtkPointsList, int& positive, int& negative) const;

			static void sortCoplanarPointsByAngle(std::vector<Point>& points, Point& centroid, Point& fixPlaneNormalIn, bool clockwise = false);

			static Plane getBestPlane(const std::vector<PointTypeITK>& pPoints, bool& result);

			static Plane getBestPlane(const std::vector<Point>& pPoints, bool& result);

			static Plane getBestPlane(const std::list<Point>& pPoints, bool& result);

			static Plane getBestPlaneOutliers(const std::vector<Point>& pPoints, std::vector<Point>& pOutLiers, bool& result, double pPercent = 0.25);

			static Plane getBestPlane(const std::vector<Point>& pPoints, const Point& pParallelVector, bool& result);

			Point getNormalVector() const;

			double getBias() const;

			Point getPoint() const;

			double eval(const Point& a) const;

			double eval(const double a[3]) const;

			void reverse();

			void reverseByNormal(const Point& checkNormal);

			void reverseByPoint(const Point& pPoint, bool sameDirection = true);

			void movePlane(const Point& pPoint);

			void movePlaneOnNormal(double distance);

			bool getIsInit() const;

			void show() const;

			void deletePlane();

		protected:
			friend class FemurImplant;
			friend class TibiaImplant;
			friend class FemurImplantMatch;
			friend class TibiaImplantMatch;
		private:
			Point normalVector;
			Point mPoint;
			double bias;
			bool isInit;
			void fixNormalVector(const Point& newNormalVector);
			void setPoint(const Point& pPoint);
		};
	}
}

#endif