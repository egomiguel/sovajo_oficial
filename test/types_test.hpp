
#ifndef TYPES_TEST_H
#define TYPES_TEST_H

#include "itkImage.h"

using SegmentPixelTypeTest = uint8_t;
using PixelTypeTest = int16_t;
using ImageTypeTest = itk::Image<PixelTypeTest, 3>;
using SegmentImageTypeTest = itk::Image<SegmentPixelTypeTest, 3>;

using PointTypeITK = itk::Point<double, 3>;

#endif