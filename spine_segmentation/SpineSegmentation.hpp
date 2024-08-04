#ifndef SPINE_SEGMENTATION_H
#define SPINE_SEGMENTATION_H

#include <string>
#include "Types.hpp"
#include "spine_segmentation_export.h"

namespace SPINE
{
	namespace SEGMENTATION
	{
		class SPINE_SEGMENTATION_EXPORT SpineSegmentation
		{
		public:

			struct Plane
			{
				itk::Point<double, 3> normal;
				itk::Point<double, 3> center;
			};

			SpineSegmentation(ImageType::Pointer pSpine);

			std::vector<SpineSegmentation::Plane> getIntervertebralPlanes(const std::vector<ImageType::PointType>& centerPhysicalPoints, const ImageType::RegionType& region) const;

		private:
			double getPixel(const ImageType::Pointer image, int x, int y, int z) const;

			ImageType::Pointer mSpine;

			double calculateMean(const std::vector<double>& values) const;

			double calculateCovariance(const std::vector<double>& x, const std::vector<double>& y, double meanX, double meanY) const;

			double calculateVariance(const std::vector<double>& values, double mean) const;

			double calculateCorrelation(const std::vector<double>& values, int lag) const;

			int getLowerValueAsPos(const std::vector<double>& values, int pBegin, int pEnd) const;

			int getGreaterValueAsPos(const std::vector<double>& values, int pBegin, int pEnd) const;

			cv::Point3d getFixNormal(const cv::Point3d& V1, const cv::Point3d& V2, bool fixFirst = true) const;

		};
	}
}

#endif