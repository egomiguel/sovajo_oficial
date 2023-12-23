#ifndef IMPLANT_UKA_BALANCE_H
#define IMPLANT_UKA_BALANCE_H

#include "Knee.hpp"
#include "LegAngle.hpp"
#include "itkImage.h"
#include "Types.hpp"
#include <itkRigid3DTransform.h>
#include "uka_implants_export.h"
#include "FemurImplant.hpp"
#include "TibiaImplant.hpp"
#include "vtkPlane.h"

namespace UKA
{
	namespace IMPLANTS
	{

		struct UKA_IMPLANTS_EXPORT BalanceInfo
		{
			double distanceLateral;
			double distanceMedial;
			BalanceInfo(double distanceLateral, double distanceMedial)
			{
				this->distanceLateral = distanceLateral;
				this->distanceMedial = distanceMedial;
			}
		};

		struct UKA_IMPLANTS_EXPORT AnglesInfo
		{
			double varus_angle;
			double flexion_angle;
			double rotation_angle;
			AnglesInfo(double flexion_angle, double varus_angle, double rotation_angle)
			{
				this->flexion_angle = flexion_angle;
				this->varus_angle = varus_angle;
				this->rotation_angle = rotation_angle;
			}
		};

		class UKA_IMPLANTS_EXPORT Balance
		{
		public:
			Balance();

			~Balance();

			void init(const Knee& pKnee, const FemurImplant& pFemurImplant, const TibiaImplant& pTibiaImplant, const itk::Rigid3DTransform<>::Pointer pImplantToBoneFemurTransform, const itk::Rigid3DTransform<>::Pointer pImplantToBoneTibiaTransform);

			void setTransformFemurCtToMarker(const itk::Rigid3DTransform<>::Pointer transform);

			void setTransformFemurMarkerToCamera(const itk::Rigid3DTransform<>::Pointer transform);

			void setTransformTibiaCtToMarker(const itk::Rigid3DTransform<>::Pointer transform);

			void setTransformTibiaMarkerToCamera(const itk::Rigid3DTransform<>::Pointer transform);

			void setTransformKneeCapCtToMarker(const itk::Rigid3DTransform<>::Pointer transform);

			void setTransformKneeCapMarkerToCamera(const itk::Rigid3DTransform<>::Pointer transform);

			BalanceInfo distanceByAngleBeforeResectionBone() const;

			//BalanceInfo distanceByAngleAfterResectionBone(double toolSize = 0) const;

			//BalanceInfo distanceByAngleAfterResectionBone(const Plane& pTibia, const vtkSmartPointer<vtkPolyData> pFemurPoly, double toolSize = 0) const;

			double distanceByAngleFemurImplantToTibiaPlane() const;

			Plane ComputeNewPlaneTibia(const itk::Rigid3DTransform<>::Pointer pImplantToBoneTibiaTransform);

			//vtkSmartPointer<vtkPolyData> ComputeNewPolyFemur(const itk::Rigid3DTransform<>::Pointer pImplantToBoneFemurTransform);

			double distanceByAngleWithImplant(const PointTypeITK& implantPlateau, bool useImplantThickness = false) const;

			AnglesInfo anglesByMotion() const;

			//itk::Rigid3DTransform<>::Pointer getNewImplantToBoneFemurTransformAxial(const std::vector<PointTypeITK>& axialPointsCT);

			itk::Rigid3DTransform<>::Pointer getNewImplantToBoneFemurTransformCoronal(const std::vector<PointTypeITK>& coronalPointsCT);

			itk::Rigid3DTransform<>::Pointer getNewImplantToBoneTibiaTransformAxial(const std::vector<PointTypeITK>& axialPointsCT);

			/////////////////////////////Test

			//vtkSmartPointer<vtkPolyData> getFemurPolyCut() const;


		private:

			Knee knee_;

			FemurImplant femurImplant;

			TibiaImplant tibiaImplant;

			vtkSmartPointer<vtkPolyData> femurImplantPoly;// , femurPolyCut;

			bool isInit;

			Plane PlaneTibia;

			cv::Mat femurTransformCtToMarker, tibiaTransformCtToMarker, femurTransformMarkerToCamera, tibiaTransformMarkerToCamera;

			cv::Mat mFemurTransformImplantToBone, mTibiaTransformImplantToBone;

			Point transformFemurPointToCamera(const Point& phisicalPoint) const;

			Point transformTibiaPointToCamera(const Point& phisicalPoint) const;

			Point transformFemurVectorToCamera(const Point& pLandmarkVector) const;

			Point transformTibiaVectorToCamera(const Point& pLandmarkVector) const;

			Plane transformTibiaPlane(const Plane& tibiaCT) const;

			ImplantImageType::PointType cvPointToITK(const Point& phisicalPoint) const;

			Point itkPointToCV(const ImplantImageType::PointType& phisicalPoint) const;

			cv::Mat Rigid3DTransformToCV(const itk::Rigid3DTransform<>::Pointer transform) const;

			cv::Mat Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer transform) const;

			cv::Mat Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer transform) const;

			Plane TransformImplantPlaneToBone(const Plane& plane, const itk::Rigid3DTransform<>::Pointer transform) const;

			vtkSmartPointer<vtkPolyData> TransformPolyFemurToCamera(const vtkSmartPointer<vtkPolyData> poly) const;

			vtkSmartPointer<vtkPolyData> TransformPolyTibiaToCamera(const vtkSmartPointer<vtkPolyData> poly) const;

			double closestPoint(const vtkSmartPointer<vtkPolyData> poly, double * point, Point& closest) const;

			double closestPoint(const vtkSmartPointer<vtkPolyData> poly, const Point& p, Point& closest) const;

			double GetTriangleArea(const Point& a, const Point& b, const Point& c) const;

			double GetDistancePointToSegment(const Point& S1, const Point& S2, const Point& P, Point& closestPoint) const;

			double GetClosestPointOnTriangle(const Point& a, const Point& b, const Point& c, const Point& myPoint, Point& closestPoint) const;

			cv::Point3d VtkToCv(double * P) const;
		};

	}
}

#endif