#ifndef EXCEPTION_IMPLANTS_H
#define EXCEPTION_IMPLANTS_H

#include <exception>
#include <string>
#include "tka_segmentation_export.h"

namespace TKA
{
	namespace SEGMENTATION
	{
		class TKA_SEGMENTATION_EXPORT SegmentationException : public std::exception
		{
			std::string msg_;
		public:
			SegmentationException(const std::string& msg) : msg_(msg) {}

			virtual const char* what() const noexcept override
			{
				return msg_.c_str();
			}
		};


		enum TKA_SEGMENTATION_EXPORT SegmentationExceptionCode
		{
			CAN_NOT_DEFINE_SEGMENTATION_REGION_FOR_EACH_BONE = 101,
			PLANES_FOR_MAKE_SLICES_ARE_NOT_PARALLEL,
			ALREADY_INITIALIZED_PLANE_ON_SEGMENTATION,
			POINTS_TO_DEFINE_PLANE_CAN_NOT_BE_COLLINEAR_ON_SEGMENTATION,
			PLANE_NORMAL_VECTOR_CAN_NOT_BE_ZERO_ON_SEGMENTATION,
			CAN_NOT_DETERMINE_PERPENDICULAR_PLANE_TO_PLANE_ON_SEGMENTATION
		};
	}
}

#endif
