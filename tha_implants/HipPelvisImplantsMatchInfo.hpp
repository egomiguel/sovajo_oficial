#ifndef THA_HIP_PELVIS_IMPLANT_MATCH_INFO_H
#define THA_HIP_PELVIS_IMPLANT_MATCH_INFO_H

#include "HipPelvisCupImplant.hpp"
#include "HipFemurStemImplant.hpp"
#include "HipFemurStemHeadImplant.hpp"
#include "HipPelvisCupImplantSimple.hpp"
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

			itk::Rigid3DTransform<>::Pointer setStemVersionAngle(double pStemVersionAngleDegree);

			itk::Rigid3DTransform<>::Pointer setStemVersionAngleNew(double pStemVersionAngleDegree);

			itk::Rigid3DTransform<>::Pointer setCupAngles(double pAbductionAngle, double pAnteversionAngle);

			itk::Vector< double, 3 > setCupTranslation(double pShifSuperior, double pShifLateral, double pShiftAnterior);

			itk::Vector< double, 3 > setStemHeadTranslation(double pShifSuperior, double pShifLateral, double pShiftAnterior);

			itk::Vector< double, 3 > setStemAxisShiftHip(double pShifSuperior, double pShifLateral, double pShiftAnterior);

			itk::Vector< double, 3 > setStemAxisShiftCup(double pShifSuperior, double pShifLateral, double pShiftAnterior);
			
			itk::Vector< double, 3 > matchStemToHipRotationCenter();

			itk::Vector< double, 3 > matchStemToCupRotationCenter();

			double getCupInclination() const;

			double getCupAntversion() const;

			double getCupShiftSuperior() const;

			double getCupShiftLateral() const;

			double getCupShiftAnterior() const;

			double getCupInclination(const itk::Rigid3DTransform<>::Pointer pImplantToBoneTransform) const;

			double getCupAntversion(const itk::Rigid3DTransform<>::Pointer pImplantToBoneTransform) const;

			double getCupShiftSuperior(const itk::Rigid3DTransform<>::Pointer pImplantToBoneTransform) const;

			double getCupShiftLateral(const itk::Rigid3DTransform<>::Pointer pImplantToBoneTransform) const;

			double getCupShiftAnterior(const itk::Rigid3DTransform<>::Pointer pImplantToBoneTransform) const;

			double getCupInclination(const Point& pVectorToHipCenter) const;

			double getCupAntversion(const Point& pVectorToHipCenter) const;

			double getCupShiftSuperior(const Point& pCenterOfRotation) const;

			double getCupShiftLateral(const Point& pCenterOfRotation) const;

			double getCupShiftAnterior(const Point& pCenterOfRotation) const;

			double getReamingDistance(const Point& pCenterOfRotation) const;

			double getStemAxisShiftSuperiorHip() const;

			double getStemAxisShiftLateralHip() const;

			double getStemAxisShiftAnteriorHip() const;

			double getStemAxisShiftSuperiorCup() const;

			double getStemAxisShiftLateralCup() const;

			double getStemAxisShiftAnteriorCup() const;

			double getStemHeadShiftSuperiorCup() const;

			double getStemHeadShiftLateralCup() const;

			double getStemHeadShiftAnteriorCup() const;

			double getStemVersion() const;

			double getRealStemVersion(const Point& pVectorNeckToHead) const;

			double getCombinedOffsetDistance(bool useLinerCenter = false, const Point linerCenter=Point()) const;

			double getRealCombinedOffsetDistance(const Point& pFinalCupCenter) const;

			double getHipLengthDistance(bool useLinerCenter = false, const Point linerCenter = Point()) const;

			double getRealHipLengthDistance(const Point& pFinalCupCenter) const;

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

			cv::Mat rotateStemVersionAngle(double pStemVersionAngleDegree, double alpha = 1);

		};
	}
}

#endif