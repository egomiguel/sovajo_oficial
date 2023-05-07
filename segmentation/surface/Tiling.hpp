#ifndef SEGMENTATION_TILING_H
#define SEGMENTATION_TILING_H

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "map"

namespace TKA
{
	namespace SEGMENTATION
	{

		class Tiling
		{
		private:
			int offSet;

			bool makeNormalization;

			bool closed;

			int downRot, upRot;

			vtkSmartPointer<vtkPoints> pointsDownRot, pointsUpRot, normaliceDownRot, normaliceUpRot;

			vtkSmartPointer<vtkPoints> originalDown, originalUp;

			vtkSmartPointer<vtkPoints> ExtractSortPoints(const vtkSmartPointer<vtkPolyData> polyData) const;

			std::vector<std::pair<int, int>> GetTilingPoints(const vtkSmartPointer<vtkPoints> downPoints, const vtkSmartPointer<vtkPoints> upPoints);

			vtkSmartPointer<vtkPoints> normalice(vtkSmartPointer<vtkPoints> points, double center[3], double sizeX, double sizeY);

			vtkSmartPointer<vtkTriangle> MakeCell(const std::pair<int, int>& node1, const std::pair<int, int>& node2, int offSetDown, int offSetUp);

			double GetDistance(const double p1[3], const double p2[3]) const;

			vtkSmartPointer<vtkPoints> RotatePoints(vtkSmartPointer<vtkPoints> points, vtkIdType id);

			std::pair<int, int> GetNearestPoints(const vtkSmartPointer<vtkPoints> pointsDown, const vtkSmartPointer<vtkPoints> pointsUp);

			void UpdatePath(std::vector<std::pair<int, int>>& path);

			vtkSmartPointer<vtkPolyData> PointsToContour(const vtkSmartPointer<vtkPoints> points);

			vtkSmartPointer<vtkPolyData> TriangulateContour(const vtkSmartPointer<vtkPolyData> contour);

			vtkSmartPointer<vtkTriangle> CreateCells(vtkIdType p1, vtkIdType p2, vtkIdType p3);

		public:
			Tiling(const vtkSmartPointer<vtkPolyData> sliceDown, const vtkSmartPointer<vtkPolyData> sliceUp, int offSet = 0, bool closed = true, bool makeNormalization = false);

			Tiling(const vtkSmartPointer<vtkPoints> pointsDown, const vtkSmartPointer<vtkPolyData> sliceUp, int offSet = 0, bool closed = true, bool makeNormalization = false);

			Tiling(const vtkSmartPointer<vtkPoints> pointsDown, const vtkSmartPointer<vtkPoints> pointsUp, int offSet = 0, bool closed = false, bool makeNormalization = false);

			vtkSmartPointer<vtkPolyData> MakeTiling();

			void MakeTiling(vtkSmartPointer<vtkPoints>& points, vtkSmartPointer<vtkCellArray>& cells, bool addFixPoints = true);

			void GetTriangulateCellsDown(vtkSmartPointer<vtkCellArray>& cells, int currentOffset);

			void GetTriangulateCellsUp(vtkSmartPointer<vtkCellArray>& cells, int currentOffset);

			vtkSmartPointer<vtkPoints> GetOriginalUp() const;

			vtkSmartPointer<vtkPoints> GetOriginalDown() const;
		};
	}
}

#endif