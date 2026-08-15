#ifndef TKA_REGISTRATION_NAV_FEMUR_H
#define TKA_REGISTRATION_NAV_FEMUR_H

#include "Registration.hpp"
#include "tka_registration_nav_export.h"

namespace TKA_NAV
{
	namespace REGISTRATION
	{
		class TKA_REGISTRATION_NAV_EXPORT FemurRegistration : public Registration
		{
		public:
			FemurRegistration(const vtkSmartPointer<vtkPolyData> img, const PointTypeITK& pHipCenterCT, const PointTypeITK& pKneeCenterCT, const PointTypeITK& pMedialEpicondyleCT);

			~FemurRegistration();

			bool MakeRegistration(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pHipCamera, const PointTypeITK& pKneeCenterCamera, const PointTypeITK& pMedialEpicondyleCamera, bool useRandomAlignment = false);

		private:
			PointTypeITK hipCenterCT;
			PointTypeITK kneeCenterCT;
			PointTypeITK medialEpicondyleCT;
		};
	}
}

#endif
