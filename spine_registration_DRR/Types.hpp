#ifndef SPINE_REGISTRATION_TYPES_H
#define SPINE_REGISTRATION_TYPES_H

#include "itkImage.h"

namespace SPINE
{
	namespace REGISTRATION_DRR
	{
		/*
		using PixelType3D = short;
		using RegistrationDRRPixelType = float;

		using RegistrationDRRImageType3D = itk::Image<PixelType3D, 3>;
		using RegistrationDRRInternalPixelType = itk::Image<RegistrationDRRPixelType, 3>;
		*/

		using SegmentPixelType = uint8_t;
		using PixelType = int16_t;
		using ImageType3D = itk::Image<PixelType, 3>;
		using ImageType2D = itk::Image<PixelType, 2>;
		using SegmentImageType = itk::Image<SegmentPixelType, 3>;
		using PointTypeITK = itk::Point<double, 3>;

		const PixelType ImageTypeMax = itk::NumericTraits<ImageType2D::PixelType>::max();
		const PixelType ImageTypeMin = itk::NumericTraits<ImageType3D::PixelType>::min();

	}
}

#endif