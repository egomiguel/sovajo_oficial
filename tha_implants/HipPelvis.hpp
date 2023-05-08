#ifndef HIP_PELVIS_H
#define HIP_PELVIS_H

#include "Plane.hpp"
#include "Utils.hpp"
#include <vector>
#include "Types.hpp"
#include "tha_implants_export.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

namespace THA
{
	namespace IMPLANTS
	{
		class THA_IMPLANTS_EXPORT HipPelvis
		{
		public:
			HipPelvis();
			void init(const Point& pLeftASIS, const Point& pRightASIS, const Point& pLeftPubicTubercle, const Point& pRightPubicTubercle, const Point& pLeftLesserTrochanter, const Point& pRightLesserTrochanter, const vtkSmartPointer<vtkPolyData>& pPelvis, PelvisSide pSide);

			Point getMidASIS() const;
			Point getPelvisVectorASIS() const;
			Point getPelvisVectorLateralASIS() const;
			Point getPelvisVectorAP() const;
			Point getPelvisVectorInfSup() const;

			Point getFemurVectorLatMed(const Point& pCenterOfRotation) const;
			Point getFemurVectorInfSup() const;

			double getHipLengthDistance() const;
			double getCombinedOffsetDistance() const;

			std::pair<Point, Point> getAbductionAnteversionVectorsZX(const Point& pCenterOfRotation, double pAbductionAngle, double pAnteversionAngle) const;
			Point getPubicJoin() const;

			vtkSmartPointer<vtkPolyData> getPelvisVTK() const;
			static Point getNativeCenterOfRotation(const std::vector<Point>& pPoints);

			///////////////////////////////Test

			vtkSmartPointer<vtkPolyData> getPelvisFemurCutTest();

		private:
			Point mLeftLesserTrochanter, mRightLesserTrochanter;
			Point mLeftASIS, mRightASIS;
			Point mLeftPubicTubercle, mRightPubicTubercle;
			Point mPubicJoin;
			Point mFemurAxisVector, mFemurPointInsideCenter;
			Plane mPlaneAPP;
			PelvisSide mSide;
			vtkSmartPointer<vtkPolyData> mPelvis;
			bool isInit;

			void extractFemurAxisVector();
		};
	}
}

#endif