#ifndef THA_HIP_FEMUR_STEM_IMPLANT_H
#define THA_HIP_FEMUR_STEM_IMPLANT_H

#include "tha_implants_export.h"
#include "Plane.hpp"

namespace THA
{
	namespace IMPLANTS
	{
		class THA_IMPLANTS_EXPORT HipFemurStemImplant
		{
		public:
			HipFemurStemImplant();
			void init(const Point& pTopPoint, const Point& pBasePoint, const Point& pHeadCenter, const Point& pNeckCenter);
			Point getVectorInfSup() const;
			Point getVectorLatMed() const;
			Point getVectorNeckToHead() const;
			Point getHeadCenter() const;
			Point getCanalAxisTopPoint() const;

		private:
			Point mTopPoint, mBasePoint, mHeadCenter, mNeckCenter;
			bool isInit;
		};

	}
}


#endif