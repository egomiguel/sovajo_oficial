#ifndef FEMUR_IMPLANT_THREE_PLANE_UKA_H
#define FEMUR_IMPLANT_THREE_PLANE_UKA_H

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include <itkRigid3DTransform.h>
#include "Plane.hpp"
#include "FemurImplant.hpp"
#include "uka_implants_export.h"
#include "Utils.hpp"

namespace UKA
{
	namespace IMPLANTS
	{
		class UKA_IMPLANTS_EXPORT FemurImplantThreePlane: public FemurImplant
		{
		public:
			FemurImplantThreePlane();

			virtual ~FemurImplantThreePlane() {};

			FemurImplantThreePlane(const FemurImplantThreePlane& pImplant);

			void init(const Plane& pPosterior, const Plane& pCenter, const Plane& pAnterior, const Point& pRodTopPoint, 
					  const Point& pRodBaseExtremeSide1, const Point& pRodBaseExtremeSide2,
					  const vtkSmartPointer<vtkPolyData> implantModel, const FemurImplantInfo& pImplantInfo);

			Plane getPosterior() const override;

			Plane getMidPlane() const;

			Plane getDistalPlane() const;

			double getWidthSize() const;

			Point getRodTopPoint() const;

			Point getRodTopPointProjectedOnBase() const;

			Point getRodTopPointProjectedOnBaseExterior() const;

			Plane getAnteriorPlane() const;

			Plane getCenterPlane() const;

			/*std::vector<Point> getSortPointsSide1() const;

			std::vector<Point> getSortPointsSide2() const;

			std::vector<Point> getAllSidePointsInOrder(double pOffset = 0) const;

			std::vector<Point> getAllSidePointsInOrder(cv::Mat& pRotation, cv::Mat& pTranslation, double pOffset = 0) const;

			Plane getBestPlaneToCurvePoints() const;*/

			Point getDirectVectorFemurAxis() const;

			Point getDirectVectorTEA() const;

			Point getDirectVectorAP() const;

			Point getRotationPoint() const;

		private:
			Plane mDistal, mAnterior, mCenter;
			Point mRodBasePoint1, mRodBasePoint2, mRodTopPoint;//, mRodTopPointProjectedOnBase;
			Point mVectorAP, mVectorTEA, mVectorForceLine;
		};

	}
}

#endif