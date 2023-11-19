#ifndef TIBIA_IMPLANT_PKA_H
#define TIBIA_IMPLANT_PKA_H

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include <itkRigid3DTransform.h>
#include "Plane.hpp"
#include "pka_implants_export.h"
#include "Utils.hpp"

namespace PKA
{
	namespace IMPLANTS
	{
		class PKA_IMPLANTS_EXPORT TibiaImplant
		{
		public:
			TibiaImplant();

			TibiaImplant(const TibiaImplant& pImplant);

			void init(const Point& apLinePclPoint, const Point& apLineTuberPoint, const Point& sidePoint, const Point& exteriorPoint, const TibiaImplantInfo& pImplantInfo);

			Plane getTibiaPlane() const;

			Point getTibiaNormalVector() const;

			Point getTibiaVectorAP() const;

			Point getTibiaVectorTEA() const;

			cv::Mat getTibiaKneeCenter() const;

			Point getCentralPoint() const;

			Point getExteriorPoint() const;

			Point getPointPCL() const;

			Point getPointTuber() const;

			Point getPointSide() const;

			TibiaImplantInfo getImplantInfo() const;

		private:
			Plane tibiaPlane;
			Point centralPoint, exteriorPoint, pclPoint, tuberPoint, sidePoint;
			TibiaImplantInfo mImplantInfo;
			bool isInit;
			//void fixNormalVectorTibia(const Point& fixPoint, const Point& referencePoint);
		};
	}
}

#endif