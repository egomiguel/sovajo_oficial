#ifndef THA_HIP_FEMUR_OPPSIDE_H
#define THA_HIP_FEMUR_OPPSIDE_H

#include "Types.hpp"
#include "Point.hpp"
#include "tha_implants_export.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

namespace THA
{
	namespace IMPLANTS
	{
		class THA_IMPLANTS_EXPORT HipFemurOppside
		{
		public:
			HipFemurOppside();
			void init(const Point& headCenter, const Point& canalCenter, const Point& lesserTrochanter,
				 const Point& femurKneeCenter);

			void init(const Point& headCenter, const Point& canalDistalCenter, const Point& canalProximalCenter, 
				const Point& lesserTrochanter, const Point& femurKneeCenter); 

			Point getHeadCenter() const;
			Point getLesserTrochanter() const;
			Point getCanalAxisVectorInfSup() const;
			Point getCanalAxisPoint() const;
			Point getKneeCenter() const;

		private:
			Point mHeadCenter;
			Point mLesserTrochanter;
			Point mCanalAxisVectorInfSup;
			Point mCanalAxisPoint;
			Point mKneeCenter;
			bool isInit;
		};
	}
}

#endif