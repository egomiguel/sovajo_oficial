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

			static itk::Rigid3DTransform<double>::Pointer ImageRegistration2D3D(const ImageType3D::Pointer& pFixedImage, const ImageType3D::Pointer& pMovingImage, ImageType3D::Pointer& pMovingImageOutput, int pDefaultPixel = 0);
			static SegmentImageType::Pointer GetDigitallyReconstructedRadiograph(const ImageType3D::Pointer& pImageCT, float pRotationX, float pRotationY, float pRotationZ, int pOutputSize, bool pVerbose = false);

		private:
			template<typename ImageTypeInput, typename ImageTypeOutput>
			static typename ImageTypeOutput::Pointer castImage(const typename ImageTypeInput::Pointer input);
		};
	}
}

#endif