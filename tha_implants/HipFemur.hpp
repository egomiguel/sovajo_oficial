#ifndef THA_HIP_FEMUR_H
#define THA_HIP_FEMUR_H

#include "Plane.hpp"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

namespace THA
{
	namespace IMPLANTS
	{
		class THA_IMPLANTS_EXPORT HipFemur
		{
		public:
			HipFemur();
			void init(const Point& headCenter, const Point& neckCenter, const Point& canalDistalCenter, const Point& canalProximalCenter, const Point& lesserTrochanter,
				const Point& medialEpicondyle, const Point& lateralEpicondyle, const Point& femurKneeCenter, 
				const vtkSmartPointer<vtkPolyData>& femurPoly);
			Point getNeckCenter() const;
			vtkSmartPointer<vtkPolyData> getFemur() const;
			Point getHeadCenter() const;
			Point getLesserTrochanter() const;
			Point getMedialEpicondyle() const;
			Point getLateralEpicondyle() const;
			Point getCanalAxisVectorInfSup() const;
			Point getVectorLatMed() const;
			Point getCanalAxisPoint() const;
			Point getProximalCanalAxisPoint() const;
			Point getNeckAxisVectorToHead() const;
			Point getKneeCenter() const;
			double getFemurVersion(const Point& pNeckAxisVectorToHead, const PelvisSide& pOperationSide) const;

		private:
			vtkSmartPointer<vtkPolyData> mFemur;
			Point mHeadCenter;
			Point mNeckCenter;
			//Point mGreaterTrochanter;
			Point mLesserTrochanter;
			Point mMedialEpicondyle;
			Point mLateralEpicondyle;
			Point mCanalAxisVectorInfSup;
			Point mNeckAxisVectorToHead;
			Point mCanalAxisPoint, mCanalAxisProximalPoint, mCanalAxisDistalPoint;
			Point mKneeCenter;
			void getNeckAxisVector();
			bool isInit;
		};
	}
}

#endif