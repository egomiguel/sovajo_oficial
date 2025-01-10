#ifndef THA_HIP_FEMUR_STEM_HEAD_IMPLANT_H
#define THA_HIP_FEMUR_STEM_HEAD_IMPLANT_H

#include "tha_implants_export.h"
#include "Plane.hpp"

namespace THA
{
	namespace IMPLANTS
	{
		class THA_IMPLANTS_EXPORT HipFemurStemHeadImplant
		{
		public:
			HipFemurStemHeadImplant();
			void init(const Point& pHeadBasePoint1, const Point& pHeadBasePoint2, const Point& pHeadBasePoint3, 
					  const Point& pHeadInsideCenterTopPoint, const std::vector<Point>& pExternalPoints);
			void init(const Point& pHeadBasePoint1, const Point& pHeadBasePoint2, const Point& pHeadBasePoint3,
				const Point& pHeadInsideCenterTopPoint, const Point& pSphereCenter, double pRadius);
			Point getVectorInfSup() const;
			Point getInsideCenterTopPoint() const;
			Point getCenterOfSphere() const;
			double getRadius() const;

		private:
			Point mHeadInfSupVector;
			Point mHeadInsideCenterTopPoint;
			Point mSphereCenter;
			double mRadius;
			bool isInit;
		};

	}
}


#endif