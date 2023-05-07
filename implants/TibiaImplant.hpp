#ifndef TIBIA_IMPLANT_H
#define TIBIA_IMPLANT_H

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include <itkRigid3DTransform.h>
#include "Plane.hpp"
#include "tka_implants_export.h"
#include "Utils.hpp"

namespace TKA
{
	namespace IMPLANTS
	{
		class TKA_IMPLANTS_EXPORT TibiaImplant
		{
		public:
			TibiaImplant();

			TibiaImplant(const TibiaImplant& pImplant);

			void init(const Point& pclPoint1, const Point& pclPoint2, const Point& frontPoint, const Point& exteriorPoint, const TibiaImplantInfo& pImplantInfo);

			Plane getTibiaPlane() const;

			Point getTibiaNormalVector() const;

			Point getTibiaVectorAP() const;

			Point getTibiaVectorTEA() const;

			cv::Mat getTibiaKneeCenter() const;

			Point getCentralPoint() const;

			Point getExteriorPoint() const;

			TibiaImplantInfo getImplantInfo() const;

		private:
			Plane tibiaPlane;
			Point centralPoint, midPoint, exteriorPoint;
			TibiaImplantInfo mImplantInfo;
			bool isInit;
			//void fixNormalVectorTibia(const Point& fixPoint, const Point& referencePoint);
		};
	}
}

#endif