#ifndef EXCEPTION_HIP_H
#define EXCEPTION_HIP_H

#include "tka_hip_export.h"

namespace TKA
{
	namespace HIP
	{

		enum TKA_HIP_EXPORT HipExceptionCode
		{
			YOU_HAVE_NOT_GENERATED_ENOUGH_DATA = 301,
			ELLIPSES_HAVE_NOT_BEEN_INITIALIZED,
			SPHERE_HAVE_NOT_BEEN_INITIALIZED
		};

	}
}

#endif
