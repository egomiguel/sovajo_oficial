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

namespace THA
{
	namespace IMPLANTS
	{
		class THA_IMPLANTS_EXPORT HipPelvis
		{
		public:
			HipPelvis();
			void init(const Point& pLeftASIS, const Point& pRightASIS, const Point& pLeftPubicTubercle, const Point& pRightPubicTubercle, const vtkSmartPointer<vtkPolyData>& pPelvis,
					  const HipFemur& pFemur, const HipFemurOppside& pFemurOppside, PelvisSide pSide);

			//Point getMidASIS() const;
			Point getRightASIS() const;
			Point getLeftASIS() const;

			Point getPelvisVectorASIS() const;
			Point getPelvisVectorLateralASIS() const;
			Point getPelvisVectorAP() const;
			Point getPelvisVectorInfSup() const;

			Plane getCoronalPlaneAPP() const;

			double getHipLengthDistance() const; 
			double getHipLengthDistanceOppsite() const;

			double getCombinedOffsetDistance() const;
			double getCombinedOffsetDistanceOppsite() const;

			double getHipLengthDistance(const Point& pMechanicalAxis) const;
			double getCombinedOffsetDistance(const Point& pCanalAxis, const Point& pCanalAxisPoint) const;

			double getFemurVersion() const;

			double getFemurVersion(const Point& pNeckAxisVectorToHead) const;

			std::pair<Point, Point> getAbductionAnteversionVectorsZX(const Point& pCenterOfRotation, double pAbductionAngle, double pAnteversionAngle) const;
			Point getPubicJoin() const;
			PelvisSide getSide() const;

			HipFemur getFemurOperationSide() const;
			HipFemurOppside getFemurOppsite() const;

			vtkSmartPointer<vtkPolyData> getPelvisVTK() const;
			static Point getNativeCenterOfRotation(const std::vector<Point>& pPoints);

			///////////////////////////////Test

			vtkSmartPointer<vtkPolyData> getPelvisFemurCutTest();

		private:
			Point mLeftASIS, mRightASIS;
			Point mLeftPubicTubercle, mRightPubicTubercle;
			Point mPubicJoin;
			Plane mPlaneAPP;
			PelvisSide mSide;
			vtkSmartPointer<vtkPolyData> mPelvis;
			HipFemur mFemurOperationSide;
			HipFemurOppside mFemurOppsite;
			bool isInit;
		};
	}
}

#endif