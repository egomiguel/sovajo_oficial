#ifndef FEMUR_REGISTRATION_H
#define FEMUR_REGISTRATION_H

#include "Registration.hpp"
#include "pka_registration_export.h"

namespace PKA
{
	namespace REGISTRATION
	{
		class PKA_REGISTRATION_EXPORT FemurRegistration : public Registration
		{
		public:
			FemurRegistration(const vtkSmartPointer<vtkPolyData> img, const PointTypeITK& pHipCenterCT, const PointTypeITK& pKneeCenterCT, const PointTypeITK& pEpicondyleCT);

			~FemurRegistration();

			bool MakeRegistration(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pHipCamera, const PointTypeITK& pKneeCenterCamera, const PointTypeITK& pEpicondyleCamera, bool useRandomAlignment = false);

		private:
			PointTypeITK hipCenterCT;
			PointTypeITK kneeCenterCT;
			PointTypeITK epicondyleCT;
		};
	}
}

#endif
