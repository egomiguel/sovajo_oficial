#ifndef GENERAL_REGISTRATION_H
#define GENERAL_REGISTRATION_H

#include "tha_registration_export.h"
#include "Types.hpp"
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/opencv.hpp>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <itkRigid3DTransform.h>

namespace THA
{
	namespace RIGISTRATION
	{
		class THA_REGISTRATION_EXPORT GeneralRegistrationPointsToPolydata
		{
		public:
			GeneralRegistrationPointsToPolydata(const vtkSmartPointer<vtkPolyData>& pPoly, const std::vector<PointTypeITK>& pTargetPointsOnPoly, const std::vector<PointTypeITK>& pSourceExternalPoints);

			itk::Rigid3DTransform<double>::Pointer MakeFinalAlignment(const std::vector<PointTypeITK>& pAlignmentPoints, double& pError);

			static itk::Rigid3DTransform<double>::Pointer MultiResImageRegistration(const RegistrationImageType::Pointer& pFixedImage, const RegistrationImageType::Pointer& pMovingImage, RegistrationImageType::Pointer& pMovingImageOutput, int pDefaultPixel = 0);

		private:
			cv::Mat mTransform;
			vtkSmartPointer<vtkPolyData> mPoly;

			template<typename ImageTypeInput, typename ImageTypeOutput>
			static typename ImageTypeOutput::Pointer castImage(const typename ImageTypeInput::Pointer input);
		};
	}
}

#endif