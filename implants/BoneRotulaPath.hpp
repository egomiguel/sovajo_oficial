#ifndef IMPLANT_BONE_PATELLA_PATH_H
#define IMPLANT_BONE_PATELLA_PATH_H

#include "Knee.hpp"
#include "Plane.hpp"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "CoordenateSystemFemur.hpp"
#include "Types.hpp"

namespace TKA
{
	namespace IMPLANTS
	{

		class BoneRotulaPath
		{
		private:
			std::vector<Point> mKneeGroovePath;
			//CoordenateSystemFemur * femurCoordenate;
			vtkSmartPointer<vtkPolyData> femurPoly;
		public:
			BoneRotulaPath(const Knee& knee);
			~BoneRotulaPath();
			std::vector<PointTypeITK> getKneeCapPathPoints() const;

			/////////////////////////////////Test
			vtkSmartPointer<vtkPolyData> CreateSphereTest(const Point& pPoint);
			vtkSmartPointer<vtkPolyData> JoinSagitalPointTest();
		};
	}
}

#endif