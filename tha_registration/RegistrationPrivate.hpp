#ifndef REGISTRATION_PRIVATE_H
#define REGISTRATION_PRIVATE_H


//#include <pcl/common/common_headers.h>
//#include <pcl/registration/icp.h>
//#include <pcl/kdtree/kdtree_flann.h>
//#include <pcl/point_cloud.h>
//#include <pcl/point_types.h>
#include <Eigen/Dense>
#include "Types.hpp"
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/opencv.hpp>
#include "vtkImplicitPolyDataDistance.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"

namespace THA
{
	namespace RIGISTRATION
	{
		class RegistrationPrivate
		{
		private:
			void getPointsAtDistanceRecursiveVTK(vtkIdType baseId, const cv::Point3d& pCenterPoint, double radius, const vtkSmartPointer<vtkPolyData> pPoly, std::vector<std::pair<cv::Point3d, cv::Point3d>>& result, std::set<vtkIdType>& visited);
		public:
			RegistrationPrivate();

			//pcl::PointCloud<pcl::PointXYZ>::Ptr PointsCT;

			std::vector<cv::Point3d> CvPointsCT;

			Eigen::Matrix4d mMatrix;

			double mError;

			//void alignCloud(const pcl::PointCloud<pcl::PointXYZ>::Ptr source_cloud, const pcl::PointCloud<pcl::PointXYZ>::Ptr target_cloud, double near_radius);

			//double getFitnessScore_max(const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_sc, const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_targe, bool show = false);

			//double getFitnessScore_max(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_targe);

			double getFitnessScore_max(const vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance, const std::vector<Eigen::Vector4d>& cloud_targe);

			Eigen::Vector3d itkPointToEigen(const PointTypeITK& point);

			//pcl::PointXYZ itkPointToPCL(const PointTypeITK& point);

			cv::Point3d itkPointToCV(const PointTypeITK& point);

			//PointTypeITK pclPointToItk(const pcl::PointXYZ& point);

			PointTypeITK cvPointToItk(const cv::Point3d& point);

			PointTypeITK vtkPointToItk(const double point[3]);

			Eigen::Vector3d cvPointToEigen(const cv::Point3d& point);

			cv::Mat cvPointToMat(const cv::Point3d& point);

			//void itkPointVectorToPCLPointVector(const std::vector<PointTypeITK>& pointsIn, pcl::PointCloud<pcl::PointXYZ>::Ptr pointsOut);

			void itkPointVectorToEigenPointVector4(const std::vector<PointTypeITK>& pointsIn, std::vector<Eigen::Vector4d>& pointsOut);

			void transformVectorEigen4(const std::vector<Eigen::Vector4d>& pointsIn, std::vector<Eigen::Vector4d>& pointsOut, const Eigen::Matrix4d& pMatrix);

			/*void reduceBiggerCloud(const pcl::PointCloud<pcl::PointXYZ>::Ptr biggerCloud, const pcl::PointCloud<pcl::PointXYZ>::Ptr smallerCloud,
				pcl::PointCloud<pcl::PointXYZ>::Ptr reduceBiggerCloudOut, double radius = 80.0);*/

				/*void reduceBiggerCloud(const pcl::PointCloud<pcl::PointXYZ>::Ptr biggerCloud, const pcl::PointCloud<pcl::PointXYZ>::Ptr smallerCloud,
					Eigen::MatrixXd& outPut, double radius);*/

					/*pcl::PointCloud<pcl::PointXYZ>::Ptr getNearPointsToReferencePoints(const std::vector<pcl::PointXYZ>& points, std::vector<pcl::PointXYZ>& point_list_1, std::vector<pcl::PointXYZ>& point_list_2,
						std::vector<pcl::PointXYZ>& point_list_3, std::vector<pcl::PointXYZ>& point_list_4, double radius = 10.0);*/

						/*void getMostNearPointsToReferencePoints(const std::vector<pcl::PointXYZ>& points, std::vector<pcl::PointXYZ>& point_list_1, std::vector<pcl::PointXYZ>& point_list_2,
							std::vector<pcl::PointXYZ>& point_list_3, std::vector<pcl::PointXYZ>& point_list_4, double radius = 10.0);*/

							//PointTypeITK getNearPointRandom(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, PointTypeITK& point, double radius);

							//PointTypeITK getNearPointRandom(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, const pcl::PointCloud<pcl::PointXYZ>::Ptr pointList, PointTypeITK& point);

			Eigen::Matrix4d InitRegistration(const std::vector<PointTypeITK>& source, const std::vector<PointTypeITK>& target);

			//PointTypeITK getRandomPoint(const pcl::PointCloud<pcl::PointXYZ>::Ptr pointList, const PointTypeITK& point);

			double GetNearPoint(const vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance, const PointTypeITK& point, PointTypeITK& nearPoint);

			double TransformPointDistance(const vtkSmartPointer<vtkImplicitPolyDataDistance> implicitDistance, const Eigen::Matrix4d& matrix, const PointTypeITK& point, PointTypeITK& pointOut);

			double GetITKPointDistance(const PointTypeITK& point1, const PointTypeITK& point2);

			void GetKneeCapPoints(const vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance, const cv::Point3d& a, const cv::Point3d& b, const cv::Point3d& vector, std::vector<PointTypeITK>& pPoints);

			double GetDistance(const double p1[3], const double p2[3]) const;

			double GetDistance(const cv::Point3d& p1, const cv::Point3d& p2) const;

			cv::Mat GetRotateZ(const cv::Point3d& vector);

			cv::Point3d TransformPoint(const cv::Point3d& pPoint, const cv::Mat& pTransform);

			double getACos(double pValue);

			std::vector<cv::Point3d> getPointsAtDistanceVTK(vtkIdType centerId, double radius, const vtkSmartPointer<vtkPolyData> pPoly);

			std::pair<std::vector<cv::Point3d>, std::vector<cv::Point3d>> getAxisPelvisVectors(const cv::Point3d& pAnteriorCamera, const cv::Point3d& pPosteriorCamera, const cv::Point3d& pSuperiorCamera, const cv::Point3d& pAnteriorCT, const cv::Point3d& pPosteriorCT, const cv::Point3d& pSuperiorCT);

			std::vector<cv::Point3d> getAxisHipFemoral(const PointTypeITK& pAnteriorFemoralNeck, const PointTypeITK& pAnteriorDistalTrochanter, const PointTypeITK& pLateralTrochanter);

			std::vector<cv::Point3d> getAxisHipFemoral(const cv::Point3d& pAnteriorFemoralNeck, const cv::Point3d& pAnteriorDistalTrochanter, const cv::Point3d& pLateralTrochanter);

			Eigen::Matrix4d getRigidTransform(const std::vector<cv::Point3d>& threeVectorsSource, const std::vector<cv::Point3d>& threeVectorstarget, const cv::Point3d& pointSource, const cv::Point3d& pointTarget);

			std::pair<double, cv::Point3d> getMinCircle(const cv::Point3d& a, const cv::Point3d& b, const cv::Point3d& c);

			std::pair<double, cv::Point3d> getMinCircle(const PointTypeITK& a, const PointTypeITK& b, const PointTypeITK& c);

			cv::Mat getRotateMatrix(const cv::Point3d& axis, double angle);

			cv::Mat GetRotateTransformVectors(const cv::Point3d& pFromVector, const cv::Point3d& pToVector);

			double getAngleBetweenVectors(const cv::Point3d& a, const cv::Point3d& b);
		};
	}
}


#endif