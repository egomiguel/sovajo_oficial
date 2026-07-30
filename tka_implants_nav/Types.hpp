#ifndef TKA_IMPLANTS_NAV_TYPES_H
#define TKA_IMPLANTS_NAV_TYPES_H

#include "itkImage.h"

namespace TKA_NAV
{
	namespace IMPLANTS
	{

		enum Side { LateralSide, MedialSide };
		enum PelvisSide { LEFT_SIDE, RIGHT_SIDE };
		using ImplantPixelType = uint8_t;
		using ImplantImageType = itk::Image<ImplantPixelType, 3>;
		using PointTypeITK = itk::Point<double, 3>;
	}
}

#endif