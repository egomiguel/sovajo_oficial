#ifndef SPINE_REGISTRATION_DRR_H
#define SPINE_REGISTRATION_DRR_H

#include "spine_registration_drr_export.h"
#include "Types.hpp"
#include <itkRigid3DTransform.h>

namespace SPINE
{
	namespace REGISTRATION_DRR
	{
		class SPINE_REGISTRATION_DRR_EXPORT SpineRegistrationDRR
		{
		public:
			SpineRegistrationDRR();

			static itk::Rigid3DTransform<double>::Pointer ImageRegistration2D3D(const RegistrationImageType::Pointer& pFixedImage, const RegistrationImageType::Pointer& pMovingImage, RegistrationImageType::Pointer& pMovingImageOutput, int pDefaultPixel = 0);

		private:
			template<typename ImageTypeInput, typename ImageTypeOutput>
			static typename ImageTypeOutput::Pointer castImage(const typename ImageTypeInput::Pointer input);
		};
	}
}

#endif