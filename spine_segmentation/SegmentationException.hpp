#ifndef SPINE_EXCEPTION_IMPLANTS_H
#define SPINE_EXCEPTION_IMPLANTS_H

#include <exception>
#include <string>
#include "spine_segmentation_export.h"

namespace SPINE
{
	namespace SEGMENTATION
	{
		class SPINE_SEGMENTATION_EXPORT SegmentationException : public std::exception
		{
			std::string msg_;
		public:
			SegmentationException(const std::string& msg) : msg_(msg) {}

			virtual const char* what() const noexcept override
			{
				return msg_.c_str();
			}
		};


		enum SPINE_SEGMENTATION_EXPORT SegmentationExceptionCode
		{
			PHYSICAL_POINT_TO_INDEX_OUTSIDE_OF_IMAGE,
			CAN_NOT_DETERMINE_PERPENDICULAR_LINE,
			POINTS_TO_DEFINE_PLANE_CAN_NOT_BE_COLLINEAR,
			ALREADY_INITIALIZED_PLANE,
			CAN_NOT_DETERMINE_PERPENDICULAR_PLANE_TO_PLANE,
			CAN_NOT_DETERMINE_INTERCEPTION_LINE_TO_PLANE
		};
	}
}

#endif
