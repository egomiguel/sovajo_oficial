#ifndef SEGMENTATION_TYPES_H
#define SEGMENTATION_TYPES_H

#include "itkImage.h"

namespace THA
{
	namespace SEGMENTATION
	{
		using SegmentPixelType = uint8_t;
		using PixelType = int16_t;
		using ImageType = itk::Image<PixelType, 3>;
		using ImageType2D = itk::Image<PixelType, 2>;
		using SegmentImageType = itk::Image<SegmentPixelType, 3>;
		using PointTypeITK = itk::Point<double, 3>;

		const PixelType ImageTypeMax = itk::NumericTraits<ImageType::PixelType>::max();
		const PixelType ImageTypeMin = itk::NumericTraits<ImageType::PixelType>::min();

	}
}

#endif