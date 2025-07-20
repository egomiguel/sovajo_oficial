#ifndef HIP_FEMUR_REGISTRATION_ANTERIOR_H
#define HIP_FEMUR_REGISTRATION_ANTERIOR_H

#include "Registration.hpp"
#include "HipFemurRegistration.hpp"
#include "tha_registration_export.h"

namespace THA
{
	namespace RIGISTRATION
	{
		class THA_REGISTRATION_EXPORT HipRegistrationFemurAnterior : public HipRegistrationFemur
		{
		public:
			HipRegistrationFemurAnterior(const vtkSmartPointer<vtkPolyData> pImage, const PointTypeITK& pAnteriorFemoralNeckCT, const PointTypeITK& pAnteriorDistalTrochanterCT, const PointTypeITK& pLateralTrochanterCT, RegisterSide pSide);

			std::vector<RegistrationPointsHip> getRegistrationPointAnteriorLateral(std::vector<PointTypeITK>& pVerificationPoints, double& pError) const;

			std::vector<RegistrationPointsHip> getRegistrationPointAnterior(std::vector<PointTypeITK>& pVerificationPoints, double& pError) const;

			bool RegistrationLandmarksAnterior(const PointTypeITK& pAnteriorFemoralNeckCamera, const PointTypeITK& pAnteriorDistalTrochanterCamera, const PointTypeITK& pLateralTrochanterCamera, double& error);

			bool MakeRegistrationAnterior(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pAnteriorFemoralNeckCamera, const PointTypeITK& pAnteriorDistalTrochanterCamera, const PointTypeITK& pLateralTrochanterCamera);
		
		private:
			cv::Mat getTemplateAlignment(double& pError) const;
		};
	}
}

#endif