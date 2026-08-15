#ifndef TKA_REGISTRATION_NAV_TIBIA_H
#define TKA_REGISTRATION_NAV_TIBIA_H

#include "Registration.hpp"
#include "tka_registration_nav_export.h"

namespace TKA_NAV
{
	namespace REGISTRATION
	{
		class TKA_REGISTRATION_NAV_EXPORT TibiaRegistration : public Registration
		{
		public:
			TibiaRegistration(const vtkSmartPointer<vtkPolyData> img, const PointTypeITK& pTibiaTubercleCT, const PointTypeITK& pLateralmalleolusCT, const PointTypeITK& pMedialmalleolusCT);

			~TibiaRegistration();

			bool MakeRegistration(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pTibiaTubercleCamera, const PointTypeITK& pLateralmalleolusCamera, const PointTypeITK& pMedialmalleolusCamera, bool useRandomAlignment = false);

		private:
			PointTypeITK tibiaTubercleCT;
			PointTypeITK lateralmalleolusCT;
			PointTypeITK medialmalleolusCT;
		};
	}
}

#endif