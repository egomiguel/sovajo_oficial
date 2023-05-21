#ifndef HIP_FEMUR_STEM_IMPLANT_MATCH_H
#define HIP_FEMUR_STEM_IMPLANT_MATCH_H

#include "HipFemurStemImplant.hpp"
#include "HipPelvis.hpp"
#include <itkRigid3DTransform.h>
#include "tha_implants_export.h"
#include "HipFemur.hpp"

namespace THA
{
	namespace IMPLANTS
	{
		class THA_IMPLANTS_EXPORT HipFemurStemImplantMatch
		{
		public:
			HipFemurStemImplantMatch();

			void init(const HipPelvis& pPelvis, const HipFemurStemImplant& pImplant, const Point& pHipCenterOfRotation);

			itk::Matrix< double, 3, 3 > GetRotationMatrix() const;

			itk::Vector< double, 3 > GetTranslationMatrix() const;

			double getStemVersion(const itk::Rigid3DTransform<>::Pointer pTransform, const HipFemur& pFemur) const;

			itk::Rigid3DTransform<>::Pointer getTransform(double pStemVersionAngleDegree, const HipFemur& pFemur) const;

		private:
			HipPelvis mPelvis;
			HipFemurStemImplant mImplant;
			Point mHipCenterOfRotation;
			cv::Mat rotationMatrix;
			cv::Mat translationMatrix;

			void getRigidTransform();

			bool isInit;
		};

	}
}


#endif