#ifndef SPINE_REGISTRATION_ICPLS_H
#define SPINE_REGISTRATION_ICPLS_H

#include <opencv2/calib3d/calib3d.hpp>
#include "vtkSmartPointer.h"
#include "vtkImplicitPolyDataDistance.h"
#include "vtkPolyData.h"
#include "Types.hpp"

namespace SPINE
{
	namespace SEGMENTATION
	{
		class LeastSquaresICP
		{
		private:
			struct GaussNewton
			{
				cv::Mat A, B;
				double totalError;
				double localError;
				GaussNewton(const cv::Mat& A, const cv::Mat& B, double totalError, double localError)
				{
					this->A = A;
					this->B = B;
					this->totalError = totalError;
					this->localError = localError;
				}
			};

			struct BatchResult
			{
				cv::Mat data;
				double error;

			};

			double chi2;

			double maxError;

			cv::Mat Rx(double angle);

			cv::Mat Ry(double angle);

			cv::Mat Rz(double angle);

			cv::Mat DRx(double angle);

			cv::Mat DRy(double angle);

			cv::Mat DRz(double angle);

			cv::Mat CreatePoint(double x, double y, double z);

			cv::Mat CreatePoint(cv::Point3d Point);

			cv::Mat DF_Rx(double angleX, double angleY, double angleZ, const cv::Mat& Point);

			cv::Mat DF_Ry(double angleX, double angleY, double angleZ, const cv::Mat& Point);

			cv::Mat DF_Rz(double angleX, double angleY, double angleZ, const cv::Mat& Point);

			cv::Mat Jacobian(const cv::Mat& data, const cv::Mat& Point);

			cv::Mat SquareError(const cv::Mat& data, const cv::Mat& source, const cv::Mat& target);

			GaussNewton GetSystem(const std::vector<cv::Point3d>& target, const cv::Mat& data, int posBegin, int posEnd, double lambda = 0);

			GaussNewton GetSystem(const std::vector<PointTypeITK>& target, const cv::Mat& data, int posBegin, int posEnd, double lambda = 0);

			std::vector<cv::Point3d> GetCorrespondence(const cv::Mat& target, const cv::Mat& data);

			void shuffleSource();

			std::vector<cv::Point3d> source;

		public:
			LeastSquaresICP(const std::vector<PointTypeITK>& sourcePoints);

			LeastSquaresICP(const std::vector<cv::Point3d>& sourcePoints);

			double LeastSquares(const cv::Mat& targetOnCT, cv::Mat& data, int iterations = 200);

			double LeastSquaresRandomInit(const cv::Mat& targetOnCT, cv::Mat& data, int iterations = 200);

			cv::Mat GetRotationMatrix(double angleX, double angleY, double angleZ);

			void setChi2(double pChi2);

			void setMaxError(double pMaxError);

			double getChi2() const;

			double getMaxError() const;

			static cv::Mat GetRotationAnglesXYZ(const std::vector<cv::Point3d>& threeVectorsSource, const std::vector<cv::Point3d>& threeVectorstarget, cv::Mat& data);
		};
	}
}

#endif