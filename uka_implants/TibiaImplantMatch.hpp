#ifndef TIBIA_IMPLANT_UKA_MATCH_H
#define TIBIA_IMPLANT_UKA_MATCH_H

#include <itkRigid3DTransform.h>
#include "TibiaImplant.hpp"
#include "Knee.hpp"
#include "Plane.hpp"
#include "uka_implants_export.h"

namespace UKA
{
	namespace IMPLANTS
	{

		class UKA_IMPLANTS_EXPORT TibiaImplantMatch
		{
		public:
			struct HullPoints
			{
				std::vector<PointTypeITK> implantPoints, sidePlanePoints;
			};

			TibiaImplantMatch();

			void init(const TibiaImplant& implant, const Knee& knee);

			Plane getTibiaPlane() const;

			Point transformImplantPoint(const Point& point) const;

			HullPoints GetHullPoints(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pPlateauTransformOut, 
				itk::Rigid3DTransform<>::Pointer pSideTransformOut, double distance = 1., double distancePcl = 1., 
				double distanceSide = 1., double sidePlaneWidth = 5., double closeCurveLateral = 0.8, 
				double closeCurveMedial = 0.8, int amount = 200) const;

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

			Point finalTransformPoint(const Point& pPoint, const itk::Rigid3DTransform<>::Pointer pTransform) const;

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
