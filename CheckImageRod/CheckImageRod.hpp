#ifndef CHECK_IMAGE_ROD_H
#define CHECK_IMAGE_ROD_H

#include "Types.hpp"
#include "CheckImageRod_export.h"

namespace IMAGE_ROD
{
	class CHECKIMAGEROD_EXPORT CheckImageRod
	{
	public:
		CheckImageRod(const ImageType::Pointer pKneeImage);
		bool Execute();

	private:
		ImageType::Pointer mImage;
		ImageType::Pointer cleanImage(const ImageType::Pointer image, PixelType outsidevalue, PixelType insidevalue,
			PixelType lowerThreshold, PixelType upperThreshold);

		ImageType::Pointer OtsuMultipleThresholds(const ImageType::Pointer image, short &threshold);
	};
}

#endif