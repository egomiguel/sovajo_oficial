#ifndef HIP_FEMUR_H
#define HIP_FEMUR_H

#include "Types.hpp"
#include "Point.hpp"
#include "tha_implants_export.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

namespace THA
{
	namespace IMPLANTS
	{
		class HipFemur
		{
		public:
			HipFemur();
			void init(const Point& headCenter, const Point& neck, const Point& greaterTrochanter, const Point& lesserTrochanter,
				const Point& medialEpicondyle, const Point& lateralEpicondyle, const vtkSmartPointer<vtkPolyData>& femurPoly,
				const PelvisSide& femurSide);

			vtkSmartPointer<vtkPolyData> getFemur() const;
			Point getHeadCenter() const;
			Point getNeck() const;
			Point getGreaterTrochanter() const;
			Point getLesserTrochanter() const;
			Point getMedialEpicondyle() const;
			Point getLateralEpicondyle() const;
			PelvisSide getFemurSide() const;

		private:
			vtkSmartPointer<vtkPolyData> mFemur;
			Point mHeadCenter;
			Point mNeck;
			Point mGreaterTrochanter;
			Point mLesserTrochanter;
			Point mMedialEpicondyle;
			Point mLateralEpicondyle;
			PelvisSide mFemurSide;
			bool isInit;
		};
	}
}

#endif