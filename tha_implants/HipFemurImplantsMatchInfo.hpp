#ifndef THA_HIP_FEMUR_IMPLANT_MATCH_INFO_H
#define THA_HIP_FEMUR_IMPLANT_MATCH_INFO_H

#include "HipFemurStemImplant.hpp"
#include "HipPelvis.hpp"
#include <itkRigid3DTransform.h>
#include "tha_implants_export.h"

namespace THA
{
	namespace IMPLANTS
	{
		class THA_IMPLANTS_EXPORT HipFemurImplantsMatchInfo
		{
		public:
			HipFemurImplantsMatchInfo(const HipPelvis& pPelvis, const HipFemurStemImplant& pImplantStem, const itk::Rigid3DTransform<>::Pointer pImplantToBoneStemTransform);
			
			~HipFemurImplantsMatchInfo();

			void setStemTransform(const itk::Rigid3DTransform<>::Pointer pImplantToBoneStemTransform);

			double getStemVersion() const;

			double getCombinedOffsetDistance() const;

			double getHipLengthDistance() const;

			itk::Rigid3DTransform<>::Pointer getITKStemTransform() const;

		private:

			HipFemurStemImplant mImplantStem;

			HipPelvis mPelvis;

			cv::Mat mRotationStem, mTranslationStem;

			cv::Mat Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer transform) const;

			cv::Mat Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer transform) const;

			Plane TransformPlane(const Plane& plane, const cv::Mat& rotation, const cv::Mat& translation) const;

			itk::Rigid3DTransform<>::Pointer getITKTransform(const cv::Mat& pRotation, const cv::Mat& pTranslation) const;

		};
	}
}

#endif