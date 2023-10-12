#ifndef THA_HIP_PELVIS_IMPLANT_MATCH_INFO_H
#define THA_HIP_PELVIS_IMPLANT_MATCH_INFO_H

#include "HipPelvisCupImplant.hpp"
#include "HipFemurStemImplant.hpp"
#include "HipFemurStemHeadImplant.hpp"
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
			HipPelvisImplantsMatchInfo(const HipPelvis& pPelvis, const Point& pHipCenterOfRotation, const HipPelvisCupImplant& pImplantCup, 
				const HipFemurStemImplant& pImplantStem, const HipFemurStemHeadImplant& pImplantStemHead,
				const itk::Rigid3DTransform<>::Pointer pImplantToBoneCupTransform, 
				const itk::Rigid3DTransform<>::Pointer pImplantToBoneStemTransform,
				const itk::Rigid3DTransform<>::Pointer pImplantHeadToStemTransform);
			
			~HipPelvisImplantsMatchInfo();

			void setCupTransform(const itk::Rigid3DTransform<>::Pointer pImplantToBoneCupTransform);

			void setStemTransform(const itk::Rigid3DTransform<>::Pointer pImplantToBoneStemTransform);

			void setStemHeadTransform(const itk::Rigid3DTransform<>::Pointer pImplantHeadToStemTransform);

			itk::Matrix< double, 3, 3 > setStemVersionAngle(double pStemVersionAngleDegree);

			itk::Rigid3DTransform<>::Pointer setCupAngles(double pAbductionAngle, double pAnteversionAngle);

			itk::Vector< double, 3 > setCupTranslation(double pShifSuperior, double pShifLateral, double pShiftAnterior);

			itk::Vector< double, 3 > setStemTranslation(double pShifSuperior, double pShifLateral, double pShiftAnterior);

			itk::Vector< double, 3 > matchStemToHipRotationCenter();

			double getCupInclination() const;

			double getCupAntversion() const;

			double getCupShiftSuperior() const;

			double getCupShiftLateral() const;

			double getCupShiftAnterior() const;

			double getStemShiftSuperior() const;

			double getStemShiftLateral() const;

			double getStemShiftAnterior() const;

			double getStemVersion() const;

			double getCombinedOffsetDistance() const;//Need to fix. Waiting for requirements

			double getHipLengthDistance() const;//Nedd to fix. Waiting for requirements

			double getCoverageFraction() const;

			itk::Rigid3DTransform<>::Pointer getITKStemTransform() const;

			itk::Rigid3DTransform<>::Pointer getITKCupTransform() const;

		private:

			Point mHipCenterOfRotation;

			HipPelvisCupImplant mImplantCup;

			HipFemurStemImplant mImplantStem;

			HipFemurStemHeadImplant mImplantStemHead;

			HipPelvis mPelvis;

			cv::Mat mRotationStem, mTranslationStem;

			cv::Mat mRotationStemHead, mTranslationStemHead;

			cv::Mat mRotationCup, mTranslationCup;

			cv::Mat Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer transform) const;

			cv::Mat Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer transform) const;

			Plane TransformPlane(const Plane& plane, const cv::Mat& rotation, const cv::Mat& translation) const;

		};
	}
}

#endif