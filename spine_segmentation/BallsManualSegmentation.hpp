#ifndef BALLS_MANUAL_SEGMENTATION_H
#define BALLS_MANUAL_SEGMENTATION_H

#include "Types.hpp"
#include "spine_segmentation_export.h"


namespace SPINE
{
	namespace SEGMENTATION
	{
		const int8_t BinaryBackgroundPixel = 0;
		const int8_t BinaryNoBackgroundPixel = 1;

		class SPINE_SEGMENTATION_EXPORT BallsManualSegmentation
		{
		public:

			BallsManualSegmentation(const ImageType::Pointer pImage);

			SegmentImageType::Pointer getSegmentBallsAndCentroids(short pMinthreshold, short pMaxthreshold, short pMinSize, std::vector<itk::Point<double, 3>>& pCentroids) const;

		private:
			ImageType::Pointer mImage;
			static SegmentImageType::Pointer CastImage(const ImageType::Pointer pInput);
		};
	}
}

#endif