#ifndef REGISTRATION_H
#define REGISTRATION_H

#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/opencv.hpp>
#include <itkRigid3DTransform.h>
#include "RPlane.hpp"
#include <stdint.h>
#include <vector>
#include <string>
#include "Types.hpp"
#include "pka_registration_export.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"

namespace PKA
{
	namespace REGISTRATION
	{

		class RegistrationPrivate;

		enum RegisterSide { LEFT, RIGHT };

		struct PKA_REGISTRATION_EXPORT RegistrationPointsHip
		{
			std::vector<PointTypeITK> points;
			RegistrationPointsHip(std::vector<PointTypeITK> pPoints);
			RegistrationPointsHip(std::vector<cv::Point3d> pPoints);
		};

		class PKA_REGISTRATION_EXPORT Registration
		{
		public:

			Registration(const vtkSmartPointer<vtkPolyData> img);

			virtual ~Registration();

			PointTypeITK TransformCTPointToMarkerPointITK(const PointTypeITK& point) const;

			PointTypeITK TransformMarkerPointToCtPointITK(const PointTypeITK& point) const;

			itk::Rigid3DTransform<double>::Pointer getTransformCtToMarker() const;

			itk::Rigid3DTransform<double>::Pointer getTransformMarkerToCt() const;

			double getDistanceToBone(double x, double y, double z) const;

			double getError() const;

			//std::vector<PointTypeITK> getPointsCT() const;

			static PointTypeITK makeItkPoint(double x, double y, double z);

			static itk::Rigid3DTransform<double>::Pointer GetTransformBetweenPoints(const std::vector<PointTypeITK>& source, const std::vector<PointTypeITK>& target);

			std::pair<cv::Point3d, double> getMinCircle(const std::vector<cv::Point3d>& pPoints, int amount, std::vector<cv::Point3d>& circlePoints);

			//static void ReadImage(const std::string& path, RegistrationImageType::Pointer& imgOut);

			//pcl::PointCloud<pcl::PointXYZ>::Ptr getPointsCT();

			/*static void drawCloudRGB(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud);

			static void drawCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::PointCloud<pcl::PointXYZ>::Ptr redPoints, bool applyCentroidToRed);

			static double distanceBetweenPCLPoints(pcl::PointXYZRGB& p1, pcl::PointXYZ& p2);*/

			//void getTestTransformPoints(std::vector<PointTypeITK>& points, std::vector<int>& pos, int amount = 40);
			//void showTestResult(pcl::PointCloud<pcl::PointXYZ>::Ptr a, pcl::PointCloud<pcl::PointXYZ>::Ptr b, std::vector<int>& pos, std::string msg = "");
			//PointTypeITK TransformPointTest(PointTypeITK& point, bool putRandom = false);
			//void SetTransformMatrixTest();
			//Eigen::Matrix3f rotationTest;
			//Eigen::Vector3f translationTest;

		protected:
			RegistrationPrivate* m_data = nullptr;

			bool isVTK;

			vtkSmartPointer<vtkPolyData> poly;

			vtkSmartPointer<vtkPolyData> getContour(const vtkSmartPointer<vtkPolyData> polyData, const cv::Point3d& pNormal, const cv::Point3d& pPoint);

			PointTypeITK getPointInsideSphere(const PointTypeITK& center, double radius);

			void deletePointsInsideRadius(std::list<cv::Point3d>& points, const cv::Point3d& centerPoint, double radius = 0.1);

			int getOneVerticePosition(const std::vector<cv::Point3d>& sortPoints, const RLine& perpendicularRefL1, const RLine& perpendicularL2, double factor = 1.0);

			std::vector<cv::Point3d> reduceSortPoints(const std::vector<cv::Point3d>& sortPoints, const RPlane& plane, bool isReverse = false);

			std::vector<cv::Point3d> reduceSortPointsByAngle(const std::vector<cv::Point3d>& sortPoints, const cv::Point3d& vector, const cv::Point3d& centerPoint, double minAngle, bool isReverse = false);

			void getPointAtRadius(const std::vector<cv::Point3d>& sortPoints, std::vector<cv::Point3d>& outPoints, const cv::Point3d& centerPoint, double squareRadius, int amount);

			cv::Mat GetTranslationRotation(const std::vector<PointTypeITK>& source, const std::vector<PointTypeITK>& target);

			void MakeResult(const cv::Mat& pData, double pError);

			//std::pair<cv::Point3d, double> getMinCircle(const std::vector<cv::Point3d>& pPoints, int amount, std::vector<cv::Point3d>& circlePoints);
		};

	}
}

#endif