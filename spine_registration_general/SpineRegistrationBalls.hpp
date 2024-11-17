#ifndef SPINE_BALLS_REGISTRATION_H
#define SPINE_BALLS_REGISTRATION_H

#include "spine_registration_export.h"
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/opencv.hpp>
#include <itkRigid3DTransform.h>
#include "Types.hpp"

namespace SPINE
{
	namespace REGISTRATION
	{
		class SPINE_REGISTRATION_EXPORT SpineRegistrationBalls
		{
		public:
			struct PointsFit
			{
				cv::Point3d normal, center;
			};

			SpineRegistrationBalls(const std::vector<PointTypeITK>& pTargetBallsOnCT, const std::vector<PointTypeITK>& pSourceExternalBalls);

			itk::Rigid3DTransform<double>::Pointer GetAlignment(double& pError);

		private:
			cv::Mat mTransform;
			cv::Mat mTargetOnCT;
			std::vector<PointTypeITK> mSourceExternalBalls;

			PointsFit getBestFit(const std::vector<PointTypeITK>& pPoints) const;
			static cv::Mat getRotateMatrix(const cv::Point3d& axis, double angle);
			static cv::Mat GetGeneralRotateTransformVectors(const cv::Point3d pFromVector, const cv::Point3d pToVector);
			static double getAngleBetweenVectors(const cv::Point3d& a, const cv::Point3d& b);
		};
	}
}

#endif