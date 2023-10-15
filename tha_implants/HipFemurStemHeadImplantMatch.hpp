#ifndef THA_HIP_FEMUR_STEM_HEAD_IMPLANT_MATCH_H
#define THA_HIP_FEMUR_STEM_HEAD_IMPLANT_MATCH_H

#include "HipFemurStemHeadImplant.hpp"
#include "HipFemurStemImplant.hpp"
#include <itkRigid3DTransform.h>
#include "tha_implants_export.h"

namespace THA
{
	namespace IMPLANTS
	{
		class THA_IMPLANTS_EXPORT HipFemurStemHeadImplantMatch
		{
		public:
			HipFemurStemHeadImplantMatch();

			void init(const HipFemurStemHeadImplant& pImplantHead, const HipFemurStemImplant& pImplantStem);

			itk::Matrix< double, 3, 3 > GetRotationMatrix() const;

			itk::Vector< double, 3 > GetTranslationMatrix() const;

			itk::Rigid3DTransform<>::Pointer getStemHeadTransform() const;

		private:
			HipFemurStemImplant mImplantStem;
			HipFemurStemHeadImplant mImplantHead;

			cv::Mat rotationMatrix;
			cv::Mat translationMatrix;

			void getRigidTransform();

			itk::Rigid3DTransform<>::Pointer getITKTransform(const cv::Mat& pRotation, const cv::Mat& pTranslation) const;

			bool isInit;
		};

	}
}


#endif