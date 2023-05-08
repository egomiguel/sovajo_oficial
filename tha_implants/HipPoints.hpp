#include "ImplantTools.hpp"

namespace THA
{
	namespace IMPLANTS
	{
		class HipPoints
		{
		private:
			vtkSmartPointer<vtkPolyData> mHip;
			Plane mObliqueTransverse;
			Point mTop;

			void reduceSortPointsByRadius(const std::vector<Point>& sortPoints, double radius, std::vector<cv::Point3d>& result);
		public:
			HipPoints(const vtkSmartPointer<vtkPolyData>& pHip, const Point& pTopPoint, const Plane& pObliqueTransverse);

			void GetPointOnTop(const Point& a, const Point& c, std::vector<cv::Point3d>& pResult);
			void GetTransversal(const Point& a, const Point& b, const Point& c, std::vector<cv::Point3d>& pResult);
			void GetSagitalDown(const Point& a, const Point& b, const Point& c, std::vector<cv::Point3d>& pResult);
		};
	}
}