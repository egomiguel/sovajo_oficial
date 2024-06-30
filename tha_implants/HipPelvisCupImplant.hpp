#ifndef THA_HIP_PELVIS_CUP_IMPLANT_H
#define THA_HIP_PELVIS_CUP_IMPLANT_H

#include "tha_implants_export.h"
#include "Plane.hpp"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

namespace THA
{
	namespace IMPLANTS
	{

		class THA_IMPLANTS_EXPORT HipPelvisCupImplant
		{
		public:
			HipPelvisCupImplant();
			void init(const Point& pTopPoint, const Point& pBasePoint1, const Point& pBasePoint2, const Point& pBasePoint3, const std::vector<Point>& pExternalPoints, double pHemiSphereResolution = 0.05, double pThickness = 0);
			Point getVectorX() const;
			Point getVectorZ() const;
			Point getTopPoint() const;
			Point getCenterOfRotationImplant() const;
			Plane getBasePlane() const;
			vtkSmartPointer<vtkPolyData> getHemiSphereCup() const;
			double getHemiSphereSurfaceArea() const;
			void setHemisphereResolution(double pResolution);
			double getThickness() const;
			double getInternalRadius() const;

		private:
			Point mTopPoint;
			Point mBasePoint1, mBasePoint2, mBasePoint3;
			vtkSmartPointer<vtkPolyData> mHemiSphereCup;
			double mRadius, mThickness;
			Point mCenter;
			Plane mBasePlane;
			double mHemiSphereSurfaceArea;
			bool isInit;
		};
	}
}



#endif