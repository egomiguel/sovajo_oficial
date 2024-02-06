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

		using RegistrationPixelType = int16_t;
		using RegistrationImageType = itk::Image<RegistrationPixelType, 3>;

		using PointTypeITK = itk::Point<double, 3>;
	}
}

#endif