#ifndef GENERAL_REGISTRATION_H
#define GENERAL_REGISTRATION_H

#include "general_registration_export.h"
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/opencv.hpp>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <itkRigid3DTransform.h>
#include "Types.hpp"

namespace GENERAL
{
	namespace REGISTRATION
	{
		class GENERAL_REGISTRATION_EXPORT GeneralRegistrationPointsToPolydata
		{
		public:
			GeneralRegistrationPointsToPolydata(const vtkSmartPointer<vtkPolyData>& pPoly, const std::vector<PointTypeITK>& pTargetPointsOnPoly, const std::vector<PointTypeITK>& pSourceExternalPoints);

			itk::Rigid3DTransform<double>::Pointer MakeFinalAlignment(const std::vector<PointTypeITK>& pAlignmentPoints, double& pError);

		private:
			cv::Mat mTransform;
			vtkSmartPointer<vtkPolyData> mPoly;

		};
	}
}

#endif