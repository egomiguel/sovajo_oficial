#ifndef TKA_HIP_NAV_EXCEPTION_H
#define TKA_HIP_NAV_EXCEPTION_H

#include "tka_hip_nav_export.h"

namespace TKA_NAV
{
	namespace HIP
	{

		enum TKA_HIP_NAV_EXPORT HipExceptionCode
		{
			YOU_HAVE_NOT_GENERATED_ENOUGH_DATA = 301,
			ELLIPSES_HAVE_NOT_BEEN_INITIALIZED,
			SPHERE_HAVE_NOT_BEEN_INITIALIZED
		};

	}
}

#endif
