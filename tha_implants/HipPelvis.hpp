#ifndef THA_HIP_PELVIS_H
#define THA_HIP_PELVIS_H

#include "Utils.hpp"
#include <vector>
#include "Types.hpp"
#include "tha_implants_export.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "HipFemur.hpp"
#include "HipFemurOppside.hpp"
#include "vtkImplicitPolyDataDistance.h"

namespace THA
{
	namespace IMPLANTS
	{
		class THA_IMPLANTS_EXPORT HipPelvis
		{
		public:
			HipPelvis();
			void init(const Point& pLeftASIS, const Point& pRightASIS, const Point& pLeftPubicTubercle, const Point& pRightPubicTubercle, const vtkSmartPointer<vtkPolyData>& pPelvis,
					  const HipFemur& pFemur, const HipFemurOppside& pFemurOppside, PelvisSide pSide, const Point& pHipCenterOfRotation, const Plane& pCoronalCT);

			void setCoronalTiltAngleDegree(double pTiltAngleDegree);
			double getCoronalTiltAngleDegree() const;

			void setHipCenterOfRotation(const Point& pHipCenterOfRotation);
			Point getHipCenterOfRotation() const;
			HipPelvis getHipPelvisCopyObj(double pTiltAngleDegree = 0) const;
			//Point getMidASIS() const;
			Point getRightASIS() const;
			Point getLeftASIS() const;

			Point getPelvisVectorASIS() const;
			Point getPelvisVectorLateralASIS() const;
			Point getPelvisVectorAP() const;
			Point getPelvisVectorInfSup() const;

			Plane getCoronalPlaneAPP() const;
			cv::Mat getFemurMechanicalAlignmentRotation() const;
			cv::Mat getFemurMechanicalAlignmentRotation(const Point& pCupCenter, const Point& pTranslation) const;
			cv::Mat getFemurMechanicalAlignmentRotationOppsite() const;
			cv::Mat getFemurCanalAlignmentRotation() const;
			cv::Mat getFemurCanalAlignmentRotation(const Point& pCupCenter, const Point& pTranslation) const;
			cv::Mat getFemurCanalAlignmentRotationOppsite() const;

			double getHipLengthDistance() const; 
			double getHipLengthDistanceOppsite() const;

			double getCombinedOffsetDistance() const;
			double getCombinedOffsetDistanceOppsite() const;

			double getHipLengthDistance(const Point& pHipHeadCenter, const cv::Mat& pFemurTranslation) const;
			double getCombinedOffsetDistance(const Point& pHipHeadCenter, const cv::Mat& pFemurTranslation) const;

			double getHipLengthDistance(const Point& pHipHeadCenter) const;
			double getCombinedOffsetDistance(const Point& pHipHeadCenter) const;

			double getFemurVersionRadian() const;

			double getFemurVersionDegree() const;

			double getFemurVersionRadian(const Point& pNeckAxisVectorToHead) const;

			double getFemurVersionDegree(const Point& pNeckAxisVectorToHead) const;

			double getNeckShaftAngle() const;

			std::pair<Point, Point> getAbductionAnteversionVectorsZX(const Point& pCenterOfRotation, double pAbductionAngle, double pAnteversionAngle) const;
			Point getPubicJoin() const;
			PelvisSide getSide() const;

			HipFemur getFemurOperationSide() const;
			HipFemurOppside getFemurOppsite() const;

			vtkSmartPointer<vtkPolyData> getPelvisVTK() const;

			vtkSmartPointer<vtkImplicitPolyDataDistance> getImplicitPelvisDistance() const;

			static std::pair<cv::Point3d, double> getNativeCenterOfRotation(const std::vector<Point>& pPoints);

			///////////////////////////////Test

			vtkSmartPointer<vtkPolyData> getPelvisFemurCutTest();

			double getHipLengthDistanceTest(const HipFemur& pFemurObj, const Plane& pPlaneAPP, const Point& PelvisVectorInfSup,
				const Point& pRightASIS, const Point& pLeftASIS) const;

		private:
			Point mLeftASIS, mRightASIS;
			Point mLeftPubicTubercle, mRightPubicTubercle;
			Point mPubicJoin;
			Plane mPlaneAPP;
			PelvisSide mSide;
			vtkSmartPointer<vtkPolyData> mPelvis;
			vtkSmartPointer<vtkImplicitPolyDataDistance> mImplicitPelvisDistance;
			HipFemur mFemurOperationSide;
			HipFemurOppside mFemurOppsite;
			bool isInit;
			double coronalTiltAngle;
			Point mHipCenterOfRotation;

			double getCoronalTiltAngle(const Plane& pCoronalCT) const;
		};
	}
}

#endif