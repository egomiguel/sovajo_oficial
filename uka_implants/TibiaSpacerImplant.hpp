#ifndef TIBIA_SPACER_IMPLANT_UKA_H
#define TIBIA_SPACER_IMPLANT_UKA_H

#include "Plane.hpp"
#include "uka_implants_export.h"
#include "Utils.hpp"

namespace UKA
{
	namespace IMPLANTS
	{
		class UKA_IMPLANTS_EXPORT TibiaSpacerImplant
		{
		public:
			TibiaSpacerImplant();

			TibiaSpacerImplant(const TibiaSpacerImplant& pImplant);

			void init(const Point& apLinePclPoint, const Point& apLineTuberPoint, const Point& plateauRefPointUp, const Point& exteriorPointDown);

			Plane getSpacerPlane() const;

			Point getSpacerNormalVector() const;

			Point getSpacerVectorAP() const;

			Point getSpacerVectorTEA() const;

			Point getSpacerKneeCenter() const;

		private:
			Plane spacerPlane;
			Point kneeCenterPoint, plateauRefPointUp, pclPoint, tuberPoint, plateauRefPointDown;
			bool isInit;
			//void fixNormalVectorTibia(const Point& fixPoint, const Point& referencePoint);
		};
	}
}

#endif