#ifndef FEMUR_IMPLANT_UKA_MATCH_H
#define FEMUR_IMPLANT_UKA_MATCH_H

#include <itkRigid3DTransform.h>
#include "FemurImplant.hpp"
#include "TibiaImplant.hpp"
#include "Knee.hpp"
#include "Plane.hpp"
#include "uka_implants_export.h"
#include "vtkPolyData.h"

namespace UKA
{
	namespace IMPLANTS
	{
		class UKA_IMPLANTS_EXPORT FemurImplantMatch
		{
		public:

			enum BoneAreaOnePlane {
				KOnePlanePosterior,
				KOnePlaneAnteriorAndDistalCurve
			};

			enum BoneAreaTwoPlanes {
				KTwoPlanePosterior,
				KTwoPlaneObliquePosterior,
				KTwoPlaneAnteriorAndDistalCurve
			};

			enum BoneAreaThreePlanes {
				KThreePlanePosterior,
				KThreePlaneCenter,
				KThreePlaneAnterior
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

			~FemurImplantMatch();

			void init(FemurImplant* implantFemur, const Knee& knee);

			itk::Matrix< double, 3, 3 > GetRotationMatrix() const;

			itk::Vector< double, 3 > GetTranslationMatrix() const;

			std::vector<PointTypeITK> GetHullPointsOnePlane(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut, BoneAreaOnePlane id, double distanceSide = 0, double distanceTop = 1.0, int amount = 200) const;

			std::vector<PointTypeITK> GetHullPointsThreePlanes(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut, BoneAreaThreePlanes id, double distanceSide = 0, double distanceTop = 1.0, int amount = 200) const;

			std::vector<PointTypeITK> GetHullPointsTwoPlanes(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut, BoneAreaTwoPlanes id, double distanceSide = 0, double distanceTop = 1.0, int amount = 200) const;

		private:

			FemurImplant* implantFemur;
			Knee knee;
			cv::Mat rotationMatrix;
			cv::Mat translationMatrix;
			bool isInit;
			void getRotationMatrix();
			bool getTranslationMatrix();

			Point getPointsOnPlane(const Plane& myPlane, std::vector<Point>& points) const;

			Point getPointsOnPlane(const Plane& myPlane, std::vector<Point>& allPoints, std::vector<Point>& pointsLat, std::vector<Point>& pointsMed) const;

			Plane finalTransformPlane(const Plane& plane, const itk::Rigid3DTransform<>::Pointer pTransform) const;

			void getCurveLikeU(const std::vector<Point>& points, const Point& downPoint, const Point& lateralPoint, const Point& medialPoint, const Point& topPoint, const Plane& midPlane, const Plane& currentPlane, const cv::Mat& pRotation, std::vector<Point>& vertices, double distanceSide, double distanceTop, int amount) const;

			ConvexHullFeatures getIncreaseBorder(const std::vector<Point>& points, const Point& downPoint, const Point& lateralPoint, const Point& medialPoint, const Point& topPoint, const Plane& midPlane, const Plane& currentPlane, const cv::Mat& pRotation, double distanceSideLat, double distanceSideMed, double distanceTop, double downLatCornerOut = 0.75, double downMedCornerOut = 0.75) const;

			std::vector<PointTypeITK> increaseVectorToAmount(const std::vector<Point>& points, int amount) const;

			Point movePointAtNormal(const Point& movePoint, const Point& nextPoint, const cv::Mat& rotationZ, double distance, bool changeMove = false, bool clockWise = false) const;

		};
	}
}

#endif
