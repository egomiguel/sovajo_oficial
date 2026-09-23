#ifndef FEMUR_REGISTRATION_H
#define FEMUR_REGISTRATION_H

#include "Registration.hpp"
#include "uka_registration_export.h"

namespace UKA
{
	namespace REGISTRATION
	{
		class UKA_REGISTRATION_EXPORT FemurRegistration : public Registration
		{
		public:
			FemurRegistration(const vtkSmartPointer<vtkPolyData> img, const PointTypeITK& pHipCenterCT, const PointTypeITK& pKneeCenterCT, const PointTypeITK& pEpicondyleCT, const PointTypeITK& pDistalCondyleCT);
			FemurRegistration(const vtkSmartPointer<vtkPolyData> img, const PointTypeITK& pHipCenterCT, const PointTypeITK& pKneeCenterCT, const PointTypeITK& pMedialEpicondyleCT);

			~FemurRegistration();

			bool MakeRegistration(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pHipCamera, const PointTypeITK& pKneeCenterCamera, const PointTypeITK& pEpicondyleCamera, const PointTypeITK& pDistalCondyleCamera = {}, bool useRandomAlignment = false);

		private:
			PointTypeITK hipCenterCT;
			PointTypeITK kneeCenterCT;
			PointTypeITK epicondyleCT;
			PointTypeITK distalCondyleCT;
			bool useDistalCondyleCT;
		};
	}
}

#endif
