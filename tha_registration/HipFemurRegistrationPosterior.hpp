#ifndef HIP_FEMUR_REGISTRATION_POSTERIOR_H
#define HIP_FEMUR_REGISTRATION_POSTERIOR_H

#include "Registration.hpp"
#include "HipFemurRegistration.hpp"
#include "tha_registration_export.h"

namespace THA
{
	namespace RIGISTRATION
	{
		class THA_REGISTRATION_EXPORT HipRegistrationFemurPosterior : public HipRegistrationFemur
		{
		public:
			HipRegistrationFemurPosterior(const vtkSmartPointer<vtkPolyData> pImage, const PointTypeITK& pPosteriorFemoralNeckCT, const PointTypeITK& pPosteriorDistalTrochanterCT, const PointTypeITK& pLateralTrochanterCT, RegisterSide pSide);

			std::vector<RegistrationPointsHip> getRegistrationPointPosteriorlateral(std::vector<PointTypeITK>& pVerificationPoints, double& pError) const;

			bool RegistrationLandmarksPosterior(const PointTypeITK& pPosteriorFemoralNeckCamera, const PointTypeITK& pPosteriorDistalTrochanterCamera, const PointTypeITK& pLateralTrochanterCamera, double& error);

			bool MakeRegistrationPosterior(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pPosteriorFemoralNeckCamera, const PointTypeITK& pPosteriorDistalTrochanterCamera, const PointTypeITK& pLateralTrochanterCamera);
		
		private:
			cv::Mat getTemplateAlignment(double& pError) const;
		};
	}
}

#endif