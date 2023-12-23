#ifndef TIBIA_REGISTRATION_H
#define TIBIA_REGISTRATION_H

#include "Registration.hpp"
#include "uka_registration_export.h"

namespace UKA
{
	namespace REGISTRATION
	{
		class UKA_REGISTRATION_EXPORT TibiaRegistration : public Registration
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