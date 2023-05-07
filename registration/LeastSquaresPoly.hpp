#ifndef REGISTRATION_POLYLS_H
#define REGISTRATION_POLYLS_H

#include <opencv2/calib3d/calib3d.hpp>
#include "vtkSmartPointer.h"
#include "vtkImplicitPolyDataDistance.h"
#include "vtkPolyData.h"
#include "Types.hpp"

namespace TKA
{
	namespace REGISTRATION
	{

		class LeastSquaresPoly
		{
		private:

			double maxError;

			cv::Mat Rx(double angle);

			cv::Mat Ry(double angle);

			cv::Mat Rz(double angle);

			cv::Mat DRx(double angle);

			cv::Mat DRy(double angle);

			cv::Mat DRz(double angle);

			cv::Mat CreatePoint(double x, double y, double z);

			cv::Mat CreatePoint(cv::Point3d Point);

			double DF_Rx(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation, const cv::Mat& Point);

			double DF_Ry(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation, const cv::Mat& Point);

			double DF_Rz(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation, const cv::Mat& Point);

			double DF_Translation(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation, const cv::Mat& Point);

			double getDistance(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation, const cv::Mat& Point);

			cv::Mat getGradient(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation);

			cv::Point3d ClosestPoint(const vtkSmartPointer<vtkPolyData>& surface, double point[3]);

			std::vector<cv::Point3d> source;

			vtkSmartPointer<vtkPolyData> surface;

			vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance;

			double computeFunction(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation);

			cv::Mat getGradientDiff(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation);

		public:
			LeastSquaresPoly(const std::vector<PointTypeITK>& pSourcePoints, const vtkSmartPointer<vtkPolyData>& pSurface);

			LeastSquaresPoly(const std::vector<cv::Point3d>& pSourcePoints, const vtkSmartPointer<vtkPolyData>& pSurface);

			double LeastSquares(cv::Mat& data, double learningRate = 0.5, int iterations = 100);
		};
	}
}

#endif