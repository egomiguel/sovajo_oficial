#ifndef FEMUR_IMPLANT_H
#define FEMUR_IMPLANT_H

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include <itkRigid3DTransform.h>
#include "Plane.hpp"
#include "tka_implants_export.h"
#include "Utils.hpp"

namespace TKA
{
	namespace IMPLANTS
	{
		class TKA_IMPLANTS_EXPORT FemurImplant
		{
		public:
			FemurImplant();

			FemurImplant(const FemurImplant& pImplant);

			void init(const Plane& A, const Plane& B, const Plane& C, const Plane& D, const Plane& E,
				const Point& P3, const Point& P4, const Point& cortexPoint, const vtkSmartPointer<vtkPolyData> implantModel,
				const std::vector<PointTypeITK>& kneeCapPath, const FemurImplantInfo& pImplantInfo);

			Plane getPlaneA() const;

			Plane getPlaneB() const;

			Plane getPlaneC() const;

			Plane getPlaneD() const;

			Plane getPlaneE() const;

			Plane getMidPlane() const;

			Point getCortexPoint() const;

			cv::Mat getCortexPointMat() const;

			Point getDirectVectorFemurAxis() const;

			Point getDirectVectorTEA() const;

			Point getDirectVectorAP() const;

			Point getPointP3() const;

			Point getPointP4() const;

			Point getMidPoint() const;

			vtkSmartPointer<vtkPolyData> GetImplantModel() const;

			vtkSmartPointer<vtkPolyData> GetTransformImplantModel(const itk::Rigid3DTransform<>::Pointer transform) const;

			std::vector<Point> GetKneeCapPath() const;

			Point getMidPlaneInterceptionPoint() const;

			FemurImplantInfo getImplantInfo() const;

		private:
			Plane planeA, planeB, planeC, planeD, planeE, midPlane;
			Point P3, P4;
			Point cortexPoint;
			bool isInit;
			void fixNormalVectorC();
			vtkSmartPointer<vtkPolyData> implantModel;
			std::vector<Point> mKneeCapPath;
			void FixKneeCapPath(const std::vector<PointTypeITK>& pKneeCapPath);
			FemurImplantInfo mImplantInfo;
		};

	}
}

#endif