#ifndef FEMUR_IMPLANT_UKA_H
#define FEMUR_IMPLANT_UKA_H

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include <itkRigid3DTransform.h>
#include "Plane.hpp"
#include "uka_implants_export.h"
#include "Utils.hpp"

namespace UKA
{
	namespace IMPLANTS
	{
		class UKA_IMPLANTS_EXPORT FemurImplant
		{
		public:
			FemurImplant();

			virtual ~FemurImplant() {};

			void init(const Plane& pPosterior, const vtkSmartPointer<vtkPolyData> implantModel, 
				const FemurImplantInfo& pImplantInfo);

			virtual Plane getPosterior() const;

			vtkSmartPointer<vtkPolyData> GetImplantModel() const;

			vtkSmartPointer<vtkPolyData> GetTransformImplantModel(const itk::Rigid3DTransform<>::Pointer transform) const;

			FemurImplantInfo getImplantInfo() const;

			virtual Plane getMidPlane() const = 0;

			virtual Plane getDistalPlane() const = 0;

			virtual double getWidthSize() const = 0;

			virtual Point getRodTopPoint() const = 0;

			virtual Point getRodTopPointProjectedOnBase() const = 0;

			virtual Point getRodTopPointProjectedOnBaseExterior() const = 0;

			virtual Point getDirectVectorFemurAxis() const = 0;

			virtual Point getDirectVectorTEA() const = 0;

			virtual Point getDirectVectorAP() const = 0;

			virtual Point getRotationPoint() const = 0;

		protected:
			Plane mPosterior;
			bool isInit;
			vtkSmartPointer<vtkPolyData> mImplantModel;
			FemurImplantInfo mImplantInfo;
		};

	}
}

#endif