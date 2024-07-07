#ifndef THA_HIP_PELVIS_LINER_IMPLANT_H
#define THA_HIP_PELVIS_LINER_IMPLANT_H

#include "tha_implants_export.h"
#include "Plane.hpp"

namespace THA
{
	namespace IMPLANTS
	{

		class THA_IMPLANTS_EXPORT HipPelvisLinerImplant
		{
		public:
			HipPelvisLinerImplant();
			void init(const Point& pTopPoint, const std::vector<Point>& pInternalSpherePoints, double pThickness = 0);
			Point getTopPoint() const;
			Point getCenterOfRotationImplant() const;
			Point getCenterToTopVector() const;
			double getExternalRadius() const;

		private:
			Point mTopPoint;
			double mRadius;
			Point mCenter;
			bool isInit;
		};
	}
}



#endif