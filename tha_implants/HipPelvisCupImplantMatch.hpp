#ifndef THA_HIP_PELVIS_CUP_IMPLANT_MATCH_H
#define THA_HIP_PELVIS_CUP_IMPLANT_MATCH_H

#include "HipPelvisCupImplant.hpp"
#include "HipPelvis.hpp"
#include <itkRigid3DTransform.h>
#include "tha_implants_export.h"

namespace THA
{
	namespace IMPLANTS
	{
		class THA_IMPLANTS_EXPORT HipPelvisCupImplantMatch
		{
		public:
			HipPelvisCupImplantMatch();

			void init(const HipPelvis& pPelvis, const HipPelvisCupImplant& pImplant);

			itk::Rigid3DTransform<>::Pointer getTransform(double pAbductionAngle = 40, double pAnteversionAngle = 20, double pShifSuperior = 0, double pShifLateral = 0, double pShiftAnterior = 0, double pTiltAngle = 0) const;

			void GetRobotTransform(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut) const;

		private:
			HipPelvis mPelvis;
			HipPelvisCupImplant mImplant;
			Point mHipCenterOfRotation;
			bool isInit;
		};

	}
}


#endif