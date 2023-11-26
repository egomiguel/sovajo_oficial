#ifndef FEMUR_IMPLANT_PKA_H
#define FEMUR_IMPLANT_PKA_H

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include <itkRigid3DTransform.h>
#include "Plane.hpp"
#include "pka_implants_export.h"
#include "Utils.hpp"

namespace PKA
{
	namespace IMPLANTS
	{
		class PKA_IMPLANTS_EXPORT FemurImplant
		{
		public:
			FemurImplant();

			FemurImplant(const FemurImplant& pImplant);

			void init(const Plane& pPosterior, const Point& pRodBasePoint, const Point& pRodTopPoint,
					  const std::vector<Point>& pSortPointsSide1, const std::vector<Point>& pSortPointsSide2, 
					  const vtkSmartPointer<vtkPolyData> implantModel, const FemurImplantInfo& pImplantInfo);

			Plane getPosterior() const;

			Plane getMidPlane() const;

			Plane getDistalPlane() const;

			double getWidthSize() const;

			Point getRodTopPoint() const;

			Point getRodBasePoint() const;

			std::vector<Point> getSortPointsSide1() const;

			std::vector<Point> getSortPointsSide2() const;

			std::vector<Point> getAllSidePointsInOrder(double pOffset = 0) const;

			std::vector<Point> getAllSidePointsInOrder(cv::Mat& pRotation, cv::Mat& pTranslation, double pOffset = 0) const;

			Plane getBestPlaneToCurvePoints() const;

			Point getDirectVectorFemurAxis() const;

			Point getDirectVectorTEA() const;

			Point getDirectVectorAP() const;

			vtkSmartPointer<vtkPolyData> GetImplantModel() const;

			vtkSmartPointer<vtkPolyData> GetTransformImplantModel(const itk::Rigid3DTransform<>::Pointer transform) const;

			FemurImplantInfo getImplantInfo() const;

		private:
			Plane mPosterior, mDistal;
			Point mRodBasePoint, mRodTopPoint;
			Point mVectorAP, mVectorTEA, mVectorForceLine, mSizeMidVector;
			std::vector<Point> mSideBorder1, mSideBorder2;
			bool isInit;
			vtkSmartPointer<vtkPolyData> mImplantModel;
			FemurImplantInfo mImplantInfo;
		};

	}
}

#endif