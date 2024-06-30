#ifndef THA_HIP_PELVIS_LINER_CUP_IMPLANT_MATCH_H
#define THA_HIP_PELVIS_LINER_CUP_IMPLANT_MATCH_H

#include "HipPelvisLinerImplat.hpp"
#include "HipPelvisCupImplant.hpp"
#include "HipPelvis.hpp"
#include <itkRigid3DTransform.h>
#include "tha_implants_export.h"

namespace THA
{
	namespace IMPLANTS
	{
		class THA_IMPLANTS_EXPORT HipPelvisLinerCupImplantMatch
		{
		public:
			HipPelvisLinerCupImplantMatch();

			void init(const HipPelvisCupImplant& pCupImplant, const HipPelvisLinerImplant& pLinerImplant);

			itk::Matrix< double, 3, 3 > GetRotationMatrix() const;

			itk::Vector< double, 3 > GetTranslationMatrix() const;

			itk::Rigid3DTransform<>::Pointer getTransform() const;

			Point getLinerCenterInCup() const;

		private:
			HipPelvisCupImplant mCupImplant;
			HipPelvisLinerImplant mLinerImplant;

			cv::Mat mRotationMatrix;
			cv::Mat mTranslationMatrix;

			Point mLinerCenterInCup;

			void getRigidTransform();

			itk::Rigid3DTransform<>::Pointer getITKTransform(const cv::Mat& pRotation, const cv::Mat& pTranslation) const;

			bool isInit;
		};

	}
}


#endif