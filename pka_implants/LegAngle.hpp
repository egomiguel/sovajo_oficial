#ifndef LEG_ANGLE_H
#define LEG_ANGLE_H

#include "Plane.hpp"
#include "Types.hpp"
#include "pka_implants_export.h"

namespace PKA
{
	namespace IMPLANTS
	{
		class PKA_IMPLANTS_EXPORT LegAngle
		{

		public:
			LegAngle();

			double getAngleBetweenVectors(const Point& a, const Point& b) const;

			double getFlexionAngle(const Point& hipCenter, const Point& femurKneeCenter, const Point& lateralEpicondyle, const Point& medialEpicondyle, const Point& tibiaKneeCenter, const Point& ankleCenter, bool isRight) const;

			double getVarusAngle(const Point& femurKneeCenter, const Point& medialCondyle, const Point& femurVectorAP, const Point& femurVectorTEA, const Point& tibiaVectorTEA, const Point& tibiaVectorAP, bool isRight) const;

			double getRotationAngle(const Point& hipCenter, const Point& femurKneeCenter, const Point& femurVectorTEA, const Point& tibiaVectorTEA, const Point& tibiaVectorAP, const Point& unkle, bool isRight) const;

		private:
			Side getInclinationSide(const Point& hipCenter, const Point& kneeCenter, const Point& latEpi, const Point& medEpi, const Point& ankle) const;
			cv::Mat pointToMat(const Point& pPoint) const;
		};
	}
}

#endif