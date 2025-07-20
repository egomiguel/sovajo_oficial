#ifndef HIP_FEMUR_REGISTRATION_H
#define HIP_FEMUR_REGISTRATION_H

#include "Registration.hpp"
#include "tha_registration_export.h"

namespace THA
{
	namespace RIGISTRATION
	{
		class THA_REGISTRATION_EXPORT HipRegistrationFemur : public Registration
		{
		public:
			HipRegistrationFemur(const vtkSmartPointer<vtkPolyData> pImage, const PointTypeITK& pFemoralNeckCT, const PointTypeITK& pDistalTrochanterCT, const PointTypeITK& pLateralTrochanterCT, RegisterSide pSide);

			bool RegistrationLandmarks(const PointTypeITK& pFemoralNeckCamera, const PointTypeITK& pDistalTrochanterCamera, const PointTypeITK& pTrochanterCamera, double& error);

			bool MakeRegistration(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pFemoralNeckCamera, const PointTypeITK& pDistalTrochanterCamera, const PointTypeITK& pLateralTrochanterCamera);

		protected:
			PointTypeITK mFemoralNeck;
			PointTypeITK mDistalTrochanter;
			PointTypeITK mLateralTrochanter;
			RegisterSide mSide;
			cv::Point3d mCenter;

			cv::Mat getTemplateAlignment(double& pError) const;
		};
	}
}

#endif