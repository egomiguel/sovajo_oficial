#ifndef FEMUR_IMPLANT_PKA_MATCH_H
#define FEMUR_IMPLANT_PKA_MATCH_H

#include <itkRigid3DTransform.h>
#include "FemurImplant.hpp"
#include "Knee.hpp"
#include "Plane.hpp"
#include "pka_implants_export.h"
#include "vtkPolyData.h"

namespace PKA
{
	namespace IMPLANTS
	{
		class PKA_IMPLANTS_EXPORT FemurImplantMatch
		{
		public:

			enum BoneArea {
				KPosteriorPlane,
				KAnteriorAndDistalCurve
			};

			struct ConvexHullFeatures
			{
				Line * topLine;
				Line * downLine;
				Line * medialLine;
				Line * lateralLine;
				Point centerPoint, lateralTopPoint, medialTopPoint;
				int lateralDownPos, medialDownPos;
				std::vector<Point> convexHull;
			};

			FemurImplantMatch();

			void init(const FemurImplant& implant, const Knee& knee);

			itk::Matrix< double, 3, 3 > GetRotationMatrix() const;

			itk::Vector< double, 3 > GetTranslationMatrix() const;

			std::vector<PointTypeITK> GetHullPoints(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut, BoneArea id, double distanceSide = 0, double distanceTop = 1.0, double angleLateral = 10, double angleMedial = 15, int amount = 200) const;

		private:

			FemurImplant implant;
			Knee knee;
			cv::Mat rotationMatrix;
			cv::Mat translationMatrix;
			bool isInit;
			void getRotationMatrix();
			bool getTranslationMatrix();

			Point getPointsOnPlane(const Plane& myPlane, std::vector<Point>& points) const;

			Plane finalTransformPlane(const Plane& plane, const itk::Rigid3DTransform<>::Pointer pTransform) const;

			void getCurveLikeU(const std::vector<Point>& points, const Point& downPoint, const Point& lateralPoint, const Point& medialPoint, const Point& topPoint, const Plane& midPlane, const Plane& currentPlane, const cv::Mat& pRotation, std::vector<Point>& vertices, double distanceSide, double distanceTop, int amount) const;

			ConvexHullFeatures getIncreaseBorder(const std::vector<Point>& points, const Point& downPoint, const Point& lateralPoint, const Point& medialPoint, const Point& topPoint, const Plane& midPlane, const Plane& currentPlane, const cv::Mat& pRotation, double distanceSideLat, double distanceSideMed, double distanceTop, double downLatCornerOut = 0.75, double downMedCornerOut = 0.75) const;

			std::vector<PointTypeITK> increaseVectorToAmount(const std::vector<Point>& points, int amount) const;

			Point movePointAtNormal(const Point& movePoint, const Point& nextPoint, const cv::Mat& rotationZ, double distance, bool changeMove = false, bool clockWise = false) const;

		};
	}
}

#endif
