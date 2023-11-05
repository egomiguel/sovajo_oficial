#ifndef FEMUR_IMPLANT_MATCH_H
#define FEMUR_IMPLANT_MATCH_H

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

			FemurImplantMatch();

			void init(const FemurImplant& implant, const Knee& knee);

			itk::Matrix< double, 3, 3 > GetRotationMatrix() const;

			itk::Vector< double, 3 > GetTranslationMatrix() const;

		private:

			FemurImplant implant;
			Knee knee;
			cv::Mat rotationMatrix;
			cv::Mat translationMatrix;
			bool isInit;
			void getRotationMatrix();
			bool getTranslationMatrix();
		};
	}
}

#endif
