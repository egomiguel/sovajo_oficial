#ifndef SLICE_BORDER_H
#define SLICE_BORDER_H

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vector"
#include "Types.hpp"
#include "SPlane.hpp"
#include "tka_segmentation_export.h"

namespace TKA
{
	namespace SEGMENTATION
	{
		class SliceBorderPrivate;

		class TKA_SEGMENTATION_EXPORT SliceBorder
		{
		public:
			SliceBorder();

			~SliceBorder();

			std::vector<PointTypeITK> GetMainPoints(const vtkSmartPointer<vtkPolyData>& slice, double curvatureMaxAngle = 160.0) const;

			void AddImproveSlice(const vtkSmartPointer<vtkPolyData>& slice);

			vtkSmartPointer<vtkPolyData> GetSurface3D(bool& result, std::string& msg) const;

			vtkSmartPointer<vtkPolyData> GetSurface3D(const std::vector<vtkSmartPointer<vtkPolyData>>& polySlices, bool& result, std::string& msg) const;

			vtkSmartPointer<vtkPolyData> GetSurface3D(const std::vector<vtkSmartPointer<vtkPolyData>>& polySlices, const double normal[3], bool& result, std::string& msg) const;

			vtkSmartPointer<vtkPolyData> SmoothSurface(const vtkSmartPointer<vtkPolyData> surface);

			// axis: x = 1, y = 2, z = 3
			std::vector<vtkSmartPointer<vtkPolyData>> MakeSlicesFromPolyData(const vtkSmartPointer<vtkPolyData>& polyData, int axis, double distance = 1.0) const;

			std::vector<vtkSmartPointer<vtkPolyData>> MakeSlicesFromPolyData(const vtkSmartPointer<vtkPolyData>& polyData, const SPlane& plane1, const SPlane& plane2, double distance = 1.0) const;

			std::vector<vtkSmartPointer<vtkPolyData>> MakeSlicesFromPolyData(const vtkSmartPointer<vtkPolyData>& polyData, const SPlane& plane1, const cv::Point3d& pPoint, double distance = 1.0) const;

			/////////////////////////////////////////Test funtions

			static vtkSmartPointer<vtkPolyData> ReadPolyData(std::string name);

			static void ShowPolyData(vtkSmartPointer<vtkPolyData> poly1, vtkSmartPointer<vtkPolyData> poly2);

			static void SavePolyData(vtkSmartPointer<vtkPolyData> poly1, std::string name);

			static vtkSmartPointer<vtkPolyData> GetCenterSlice(vtkSmartPointer<vtkPolyData> polyData, int axis);

		private:
			vtkSmartPointer<vtkPolyData> RemoveBranchOnSlice(vtkSmartPointer<vtkPolyData> poly, double normal[3]);
			bool MergeTwoPolyData(const vtkSmartPointer<vtkPolyData> poly1, const vtkSmartPointer<vtkPolyData> poly2, double normal[3], vtkSmartPointer<vtkPolyData>& result) const;
			void ExtractSortLines(const vtkSmartPointer<vtkPolyData> polyData, std::list<std::pair<vtkIdType, vtkIdType>>& lines) const;
			double GetDistance(const double p1[3], const double p2[3]) const;
			SliceBorderPrivate* m_data = nullptr;
			std::vector<vtkSmartPointer<vtkPolyData>> slices;
		};
	}
}

#endif