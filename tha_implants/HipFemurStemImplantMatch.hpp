#ifndef THA_HIP_FEMUR_STEM_IMPLANT_MATCH_H
#define THA_HIP_FEMUR_STEM_IMPLANT_MATCH_H

#include "HipFemurStemImplant.hpp"
#include "HipPelvis.hpp"
#include <itkRigid3DTransform.h>
#include "tha_implants_export.h"

namespace THA
{
	namespace IMPLANTS
	{
		class THA_IMPLANTS_EXPORT HipFemurStemImplantMatch
		{
		public:
			HipFemurStemImplantMatch();

			void init(const HipPelvis& pPelvis, const HipFemurStemImplant& pImplant);

			itk::Matrix< double, 3, 3 > GetRotationMatrix() const;

			itk::Vector< double, 3 > GetTranslationMatrix() const;

			itk::Rigid3DTransform<>::Pointer getStemTransform(double pStemVersionAngleDegree) const;

			itk::Rigid3DTransform<>::Pointer getStemTransform() const;

		private:
			HipPelvis mPelvis;
			HipFemurStemImplant mImplant;
			Point mHipCenterOfRotation;
			cv::Mat rotationMatrix;
			cv::Mat translationMatrix;

			void getRigidTransform();

			itk::Rigid3DTransform<>::Pointer getITKTransform(const cv::Mat& pRotation, const cv::Mat& pTranslation) const;

			bool isInit;
		};

	}
}


#endif