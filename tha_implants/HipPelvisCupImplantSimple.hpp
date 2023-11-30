#ifndef THA_HIP_PELVIS_CUP_IMPLANT_SIMPLE_H
#define THA_HIP_PELVIS_CUP_IMPLANT_SIMPLE_H

#include "tha_implants_export.h"
#include "Plane.hpp"

namespace THA
{
	namespace IMPLANTS
	{

		class THA_IMPLANTS_EXPORT HipPelvisCupImplantSimple
		{
		public:
			HipPelvisCupImplantSimple();
			void init(const Point& pTopPoint, const Point& pBasePoint1, const Point& pBasePoint2, const Point& pBasePoint3, const std::vector<Point>& pExternalPoints);
			void init(const Point& pVectorZ, const std::vector<Point>& pExternalPoints);
			void init(const Point& pVectorZ, const Point& pCenter);
			Point getVectorZ() const;
			Point getCenterOfRotationImplant() const;

		private:
			Point mVectorZ;
			Point mCenter;
			bool isInit;
		};
	}
}



#endif