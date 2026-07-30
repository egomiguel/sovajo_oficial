#ifndef TKA_IMPLANTS_NAV_PATELLA_IMPLANT_H
#define TKA_IMPLANTS_NAV_PATELLA_IMPLANT_H

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include <itkRigid3DTransform.h>
#include "Plane.hpp"
#include "tka_implants_nav_export.h"
#include "Utils.hpp"

namespace TKA_NAV
{
	namespace IMPLANTS
	{
		class TKA_IMPLANTS_NAV_EXPORT PatellaImplant
		{
		public:
			PatellaImplant();

			PatellaImplant(const PatellaImplant& pImplant);

			void init(const Point& basePoint1, const Point& basePoint2, const Point& basePoint3, const Point& topCentralPoint);

			Plane getBasePlane() const;

			Point getNormalVector() const;

			Point getCentralPointOnBase() const;

			Point getTopPoint() const;

			double getThickness() const;

		private:
			Plane mBasePlane;
			Point mBasePoint1, mBasePoint2, mBasePoint3, mTopCentralPoint;
			bool isInit;
		};
	}
}

#endif