#ifndef TIBIA_IMPLANT_MATCH_H
#define TIBIA_IMPLANT_MATCH_H

#include <itkRigid3DTransform.h>
#include "TibiaImplant.hpp"
#include "Knee.hpp"
#include "Plane.hpp"
#include "tka_implants_export.h"

namespace TKA
{
	namespace IMPLANTS
	{

		class TKA_IMPLANTS_EXPORT TibiaImplantMatch
		{
		public:
			TibiaImplantMatch();

			void init(const TibiaImplant& implant, const Knee& knee);

			Plane getTibiaPlane() const;

			Point transformImplantPoint(const Point& point) const;

			std::vector<Point> getPointsNearImplant(double distance = 0) const;

			std::vector<PointTypeITK> GetHullPoints(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut, double distance = 1., double distancePcl = 1., int amount = 200) const;

			itk::Matrix< double, 3, 3 > GetRotationMatrix() const;

			itk::Vector< double, 3 > GetTranslationMatrix() const;

			vtkSmartPointer<vtkPolyData> GetCuttingTibia() const;

		protected:
			friend class Balance;
		private:
			TibiaImplant implant;
			Knee knee;
			cv::Mat rotationMatrix;
			cv::Mat translationMatrix;
			bool isInit;

			void makeRotationMatrix();

			void makeTranslationMatrix();

			Plane transformPlane(const Plane& plane) const;

			Plane finalTransformPlane(const Plane& plane, const itk::Rigid3DTransform<>::Pointer pTransform) const;

			Knee getKnee() const;

			bool deletePointsInsideRadius(std::list<Point>& points, const Point& centerPoint, const Point& diffPoint, Point& nearPoint, double radius = 0.2) const;

			std::vector<Point> removeOutLiers(const std::vector<Point>& points, const Point& initPoint, const Point& lastPoint) const;

			Point movePointAtNormal(const Point& movePoint, const Point& nextPoint, const cv::Mat& rotationZ, double distance, bool changeMove = false, bool clockWise = true) const;

			std::vector<PointTypeITK> increaseVectorToAmount(const std::vector<Point>& points, int amount) const;

			vtkSmartPointer<vtkPolyData> getContour(const vtkSmartPointer<vtkPolyData> poly, const Point& pNormal, const Point& pPoint) const;
		};
	}
}

#endif
