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

			void init(const Plane& pPosterior, const vtkSmartPointer<vtkPolyData> implantModel, 
				const FemurImplantInfo& pImplantInfo);

			Plane getPosterior() const;

			vtkSmartPointer<vtkPolyData> GetImplantModel() const;

			vtkSmartPointer<vtkPolyData> GetTransformImplantModel(const itk::Rigid3DTransform<>::Pointer transform) const;

			FemurImplantInfo getImplantInfo() const;

		protected:
			Plane mPosterior;
			bool isInit;
			vtkSmartPointer<vtkPolyData> mImplantModel;
			FemurImplantInfo mImplantInfo;
		};

	}
}

#endif