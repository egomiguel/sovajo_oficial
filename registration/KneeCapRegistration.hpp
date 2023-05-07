#ifndef KNEECAP_REGISTRATION_H
#define KNEECAP_REGISTRATION_H

#include "Registration.hpp"
#include "tka_registration_export.h"

namespace TKA
{
	namespace REGISTRATION
	{

		class TKA_REGISTRATION_EXPORT KneeCapRegistration : public Registration
		{
		public:
			KneeCapRegistration(const vtkSmartPointer<vtkPolyData> img, const PointTypeITK& pHipCenterCT, const PointTypeITK& pKneeCenterCT, const PointTypeITK& pLateralEpiCT, const PointTypeITK& pMedialEpiCT);

			~KneeCapRegistration();

			bool MakeRegistration(const std::vector<PointTypeITK>& pBonePoints);

			std::vector<PointTypeITK> GetRegistrationPoints() const;

		private:

			bool MakeRegistrationLS(const std::vector<PointTypeITK>& pBonePoints);

			cv::Mat GetRotateZ(const cv::Point3d& vector) const;

			cv::Point3d GetPointOnCircle(const cv::Point3d& a, const cv::Point3d& b, const cv::Point3d& center, double perCent, double radius) const;

			void GetMainPoints();

			std::vector<PointTypeITK> mPoints;

			cv::Mat transformZ;

			RPlane fitPlane;

			cv::Point3d kneeHipVector;

			cv::Point3d medLatVector;
		};
	}
}

#endif