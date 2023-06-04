#ifndef THA_HIP_PELVIS_IMPLANT_MATCH_INFO_H
#define THA_HIP_PELVIS_IMPLANT_MATCH_INFO_H

#include "HipPelvisCupImplant.hpp"
#include "HipPelvis.hpp"
#include <itkRigid3DTransform.h>
#include "tha_implants_export.h"

namespace THA
{
	namespace IMPLANTS
	{
		class THA_IMPLANTS_EXPORT HipPelvisImplantsMatchInfo
		{
		public:
			HipPelvisImplantsMatchInfo(const HipPelvis& pPelvis, const HipPelvisCupImplant& pImplantCup, const Point& pHipCenterOfRotation, const itk::Rigid3DTransform<>::Pointer pImplantToBoneCupTransform);
			
			~HipPelvisImplantsMatchInfo();

			void setCupTransform(const itk::Rigid3DTransform<>::Pointer pImplantToBoneCupTransform);

			double getCupInclination() const;

			double getCupAntversion() const;

			double getCupShiftSuperior() const;

			double getCupShiftLateral() const;

			double getCupShiftAnterior() const;

			itk::Rigid3DTransform<>::Pointer getITKCupTransform() const;

		private:

			Point mHipCenterOfRotation;

			HipPelvisCupImplant mImplantCup;

			HipPelvis mPelvis;

			cv::Mat mRotationCup, mTranslationCup;

			cv::Mat Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer transform) const;

			cv::Mat Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer transform) const;

			Plane TransformPlane(const Plane& plane, const cv::Mat& rotation, const cv::Mat& translation) const;

		};
	}
}

#endif