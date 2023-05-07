#ifndef CHECK_IMAGE_ROD_TYPES_H
#define CHECK_IMAGE_ROD_TYPES_H

#include "itkImage.h"

namespace IMAGE_ROD
{
	using SegmentPixelType = uint8_t;
	using PixelType = int16_t;
	using ImageType = itk::Image<PixelType, 3>;
	using ImageType2D = itk::Image<PixelType, 2>;
	using SegmentImageType = itk::Image<SegmentPixelType, 3>;
	using PointTypeITK = itk::Point<double, 3>;
}

#endif