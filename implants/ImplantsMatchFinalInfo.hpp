#ifndef IMPLANT_MATCH_FINAL_INFO_H
#define IMPLANT_MATCH_FINAL_INFO_H

#include "Knee.hpp"
#include <itkRigid3DTransform.h>
#include "Point.hpp"
#include "FemurImplant.hpp"
#include "TibiaImplant.hpp"
#include "tka_implants_export.h"

namespace TKA
{
	namespace IMPLANTS
	{
		class CoordenateSystemFemur;

		struct TKA_IMPLANTS_EXPORT ResectionThickness
		{
			double medial, lateral;
		};


		class TKA_IMPLANTS_EXPORT ImplantsMatchFinalInfo
		{
		public:
			ImplantsMatchFinalInfo(Knee* pKnee, const FemurImplant pFemurImplant, const TibiaImplant pTibiaImplant, const itk::Rigid3DTransform<>::Pointer pImplantToBoneFemurTransform, const itk::Rigid3DTransform<>::Pointer pImplantToBoneTibiaTransform);
			~ImplantsMatchFinalInfo();
			ResectionThickness GetTibiaResection() const;
			ResectionThickness GetFemurResectionAxial() const;
			ResectionThickness GetFemurResectionCoronal() const;

			double GetFemurVarusAngle() const;
			double GetFemurFlexionAngle() const;
			double GetFemurImplantTEAAngle() const;
			double GetFemurImplantPCAAngle() const;

			double GetTibiaVarusAngle() const;
			double GetTibiaImplantSlopeAngle() const;
			double GetTibiaImplantRotationAngle() const;

			std::pair<Point, Point> GetFemurBoneTEALine() const;
			std::pair<Point, Point> GetFemurImplantTEALine() const;
			std::vector<PointTypeITK> getBoneKneeCapPath() const;
			std::vector<PointTypeITK> getImplantKneeCapPath() const;

			/*
			   Change the medial thickness of the tibia. Updates the tibia transformation internally
			   on the object with the new translation. Return the new translation
			*/
			itk::Vector< double, 3 > SetThicknessTibiaMedial(double medialThickness);

			/*
			   Change the lateral thickness of the tibia. Updates the tibia transformation internally
			   on the object with the new translation. Return the new translation
			*/
			itk::Vector< double, 3 > SetThicknessTibiaLateral(double lateralThickness);

			/*
			   Change the medial thickness of the femur axial view. Updates the femur transformation internally
			   on the object with the new translation. Return the new translation
			*/
			itk::Vector< double, 3 > SetThicknessFemurAxialMedial(double medialThickness);

			/*
			   Change the lateral thickness of the femur axial view. Updates the femur transformation internally
			   on the object with the new translation. Return the new translation
			*/
			itk::Vector< double, 3 > SetThicknessFemurAxialLateral(double lateralThickness);

			/*
			  Change the medial thickness of the femur coronal view. Updates the femur transformation internally
			  on the object with the new translation. Return the new translation
		   */
			itk::Vector< double, 3 > SetThicknessFemurCoronalMedial(double medialThickness);

			/*
			   Change the lateral thickness of the femur coronal view. Updates the femur transformation internally
			   on the object with the new translation. Return the new translation
			*/
			itk::Vector< double, 3 > SetThicknessFemurCoronalLateral(double lateralThickness);

			/*
				Change the flexion angle of the femur implant. Updates the femur implant transformation internally
				in the object. Return the new transformation.
			*/
			const itk::Rigid3DTransform<>::Pointer setFemurFlexionAngle(double angle);

			/*
				Change the varus angle of the femur implant. Updates the femur implant transformation internally
				in the object. Return the new transformation.
			*/
			const itk::Rigid3DTransform<>::Pointer setFemurVarusAngle(double angle);

			/*
				Change the TEA rotation angle of the femur implant. Updates the femur implant transformation internally
				in the object. Return the new transformation.
			*/
			const itk::Rigid3DTransform<>::Pointer setFemurTEAAngle(double angle);

			/*
				Change the PCA angle of the femur implant. Updates the femur implant transformation internally
				in the object. Return the new transformation.
			*/
			const itk::Rigid3DTransform<>::Pointer setFemurPCAAngle(double angle);

			/*
				Change the varus angle of the tibia implant. Updates the tibia implant transformation internally
				in the object. Return the new transformation.
			*/
			const itk::Rigid3DTransform<>::Pointer setTibiaVarusAngle(double angle);

			/*
				Change the rotation angle of the tibia implant. Updates the tibia implant transformation internally
				in the object. Return the new transformation.
			*/
			const itk::Rigid3DTransform<>::Pointer setTibiaRotationAngle(double angle);

			/*
			   Change the slope angle of the tibia implant. Updates the tibia implant transformation internally
			   in the object. Return the new transformation.
		   */
			const itk::Rigid3DTransform<>::Pointer setTibiaSlopeAngle(double angle);

			/*
				Change the thickness of the tibia. Updates the tibia transformation internally
				on the object with the new translation. Return the new translation
			*/
			//itk::Vector< double, 3 > SetThicknessTibiaImplant(double lateralThickness, double medialThickness);

			/*
				Update the femur implant. Updates the transformation internally on the object for
				the new implant. Return the new transformation.
			*/
			const itk::Rigid3DTransform<>::Pointer setFemurImplant(const FemurImplant& pFemurImplant);

			/*
				Update the tibia implant. Updates the transformation internally on the object
				for the new implant. Return the new transformation.
			*/
			const itk::Rigid3DTransform<>::Pointer setTibiaImplant(const TibiaImplant& pTibiaImplant);

			/*
				Update the thickness of the femur cartilage.
			*/
			void SetCartilageFemur(double pCartilage);

			/*
				Update the thickness of the tibia cartilage.
			*/
			void SetCartilageTibia(double pCartilage);

			/*
				Update the femur transformation
			*/
			void setFemurTransform(const itk::Rigid3DTransform<>::Pointer pImplantToBoneFemurTransform);

			/*
				Update the tibia transformation
			*/
			void setTibiaTransform(const itk::Rigid3DTransform<>::Pointer pImplantToBoneTibiaTransform);

			/*
				Get the femur transformation
			*/
			itk::Rigid3DTransform<>::Pointer getITKFemurTransform() const;

			/*
				Get the tibia transformation
			*/
			itk::Rigid3DTransform<>::Pointer getITKTibiaTransform() const;

			//Knee getKnee() const;

			void test();

		private:

			Point femurAxisImplantVector, femurAPImplantVector, femurTEAImplantVector;

			Point tibiaAxisImplantVector, tibiaAPImplantVector, tibiaTEAImplantVector;

			Point femurImplantMiddlePoint, tibiaImplantMiddlePoint;

			Knee* knee;

			FemurImplant femurImplant;

			TibiaImplant tibiaImplant;

			//CoordenateSystemFemur * femurCoordenate;

			cv::Mat femurRotation, femurTranslation;

			cv::Mat tibiaRotation, tibiaTranslation;

			cv::Mat Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer transform) const;

			cv::Mat Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer transform) const;

			Plane TransformPlane(const Plane& plane, const cv::Mat& rotation, const cv::Mat& translation) const;

			void updateFemurImplantVectors();

			void updateTibiaImplantVectors();

			double femurCartilage, tibiaCartilage;

			//std::vector<PointTypeITK> boneKneeCapPath;

			std::vector<Point> implantKneeCapPath, boneKneeCapPath;
		};
	}
}

#endif