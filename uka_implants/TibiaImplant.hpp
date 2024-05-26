#ifndef TIBIA_IMPLANT_UKA_H
#define TIBIA_IMPLANT_UKA_H

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include <itkRigid3DTransform.h>
#include "Plane.hpp"
#include "uka_implants_export.h"
#include "Utils.hpp"

namespace UKA
{
	namespace IMPLANTS
	{
		class UKA_IMPLANTS_EXPORT TibiaImplant
		{
		public:
			TibiaImplant();

			TibiaImplant(const TibiaImplant& pImplant);

			void init(const Point& apLinePclPoint, const Point& apLineTuberPoint, const Point& plateauRefPointUp, const Point& exteriorPointDown, const TibiaImplantInfo& pImplantInfo);

			Plane getTibiaPlane() const;

			Point getTibiaNormalVector() const;

			Point getTibiaVectorAP() const;

			Point getTibiaVectorTEA() const;

			cv::Mat getTibiaKneeCenter() const;

			Point getCentralPoint() const;

			Point getCentralPointUp() const;

			Point getExteriorPoint() const;

			Point getPointPCL() const;

			Point getPointTuber() const;

			Point getExtremeSidePoint() const;

			Point getPlateauRefPointDown() const;

			TibiaImplantInfo getImplantInfo() const;

		private:
			Plane tibiaPlane;
			Point centralPoint, centralPointUp, exteriorPoint, pclPoint, tuberPoint, plateauRefPointDown;
			TibiaImplantInfo mImplantInfo;
			bool isInit;
			//void fixNormalVectorTibia(const Point& fixPoint, const Point& referencePoint);
		};
	}
}

#endif