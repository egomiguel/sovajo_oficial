#ifndef HIP_FEMUR_STEM_IMPLANT_H
#define HIP_FEMUR_STEM_IMPLANT_H

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
			void init(const Point& pTopPoint, const Point& pBasePoint, const Point& pHeadCenter);
			Point getVectorInfoSup() const;
			Point getVectorLatMed() const;
			Point getHeadCenter() const;

		private:
			Point mTopPoint, mBasePoint, mHeadCenter;
			bool isInit;
		};

	}
}


#endif