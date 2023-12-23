//#include <pcl/common/common_headers.h>
//#include <pcl/features/normal_3d.h>
//#include <pcl/visualization/pcl_visualizer.h>
//#include <pcl/visualization/cloud_viewer.h>
//#include <pcl/console/parse.h>
//#include <pcl/filters/uniform_sampling.h>

//#include <pcl/io/io.h>
//#include <pcl/io/pcd_io.h>
//#include <pcl/io/vtk_lib_io.h>

//#include <pcl/console/parse.h>
//#include <pcl/filters/extract_indices.h>
//#include <pcl/point_types.h>
//#include <pcl/sample_consensus/ransac.h>
//#include <pcl/sample_consensus/sac_model_plane.h>
//#include <pcl/sample_consensus/sac_model_sphere.h>
//#include <pcl/filters/voxel_grid.h>
//#include <pcl/recognition/ransac_based/trimmed_icp.h>

//#include <pcl/kdtree/kdtree_flann.h>

#include <itkBinaryMask3DMeshSource.h>
//#include <itkImageFileReader.h>

//#include <itkDirectory.h>
//#include <itkGDCMImageIO.h>
#include <itkTileImageFilter.h>
//#include <itkImageFileReader.h>
//#include <itkImageFileWriter.h>
#include <itkImage.h>
//#include <itkGDCMSeriesFileNames.h>
//#include <itkImageSeriesReader.h>
#include "vtkPlane.h"
#include "vtkCutter.h"
#include <itkRescaleIntensityImageFilter.h>
#include <iostream>
#include <thread>
#include <chrono>
#include <random>
//#include <fstream>
#include <itkMatrix.h>
#include <itkVersorRigid3DTransform.h>
#include "RegistrationException.hpp"
#include "Registration.hpp"

///////////////////////////////////////////////////////////////////////////
#include "RLine.hpp"
#include "RPlane.hpp"
#include "RegistrationPrivate.hpp"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkCleanPolyData.h"


using namespace std::chrono_literals;
using namespace UKA::REGISTRATION;

inline std::vector<PointTypeITK> PointToITKVector(const std::vector<cv::Point3d>& pData)
{
    std::vector<PointTypeITK> result;
    auto it1 = pData.begin();
    auto it2 = pData.end();

    for (; it1 != it2; ++it1)
    {
        PointTypeITK temp;
        temp[0] = (*it1).x;
        temp[1] = (*it1).y;
        temp[2] = (*it1).z;
        result.push_back(temp);
    }
    return result;
}

RegistrationPointsHip::RegistrationPointsHip(std::vector<PointTypeITK> pPoints)
{
    points = pPoints;
}

RegistrationPointsHip::RegistrationPointsHip(std::vector<cv::Point3d> pPoints)
{
    points = PointToITKVector(pPoints);
}

/*
Registration::Registration(const std::vector<RegistrationImageType::Pointer>& img)
{
    m_data = new RegistrationPrivate();
    using MeshType = itk::Mesh<double, 3>;
    auto meshSource = itk::BinaryMask3DMeshSource<RegistrationImageType, MeshType>::New();

    const auto objectValue = static_cast<RegistrationPixelType>(1);
    meshSource->SetObjectValue(objectValue);

    std::vector<MeshType::PointsContainer::ConstPointer> points;

    auto imgIt1 = img.begin();
    auto imgIt2 = img.end();
    for (; imgIt1 != imgIt2; ++imgIt1)
    {
        try
        {
            meshSource->SetInput(*imgIt1);
            meshSource->Update();
        }
        catch (itk::ExceptionObject & exp)
        {
            std::string msg = exp.what();
            throw RegistrationExceptionCode::CAN_NOT_CONVERT_ITK_IMAGE_TO_3D_MESH_SURFACE;
        }

        MeshType::ConstPointer mesh = meshSource->GetOutput();
        if (mesh->GetPoints())
        {
            if (mesh->GetPoints()->size() > 0)
            {
                points.push_back(mesh->GetPoints());
            }
            else
            {
                throw RegistrationExceptionCode::LIST_OF_IMAGE_MESH_POINTS_IS_EMPTY;
            }
        }
        else
        {
            throw RegistrationExceptionCode::CAN_NOT_GET_LIST_OF_IMAGE_MESH_POINTS;
        }
    }

    auto pointsBegin = points.begin();
    auto pointsEnd = points.end();

    itk::Point<double, 3> itkPoint;
    pcl::PointXYZ pclPoint;
    while (pointsBegin != pointsEnd)
    {
        MeshType::PointsContainer::ConstPointer bonePoints = *pointsBegin;

        auto boneIt1 = bonePoints->Begin();
        auto boneIt2 = bonePoints->End();

        for (; boneIt1 != boneIt2; ++boneIt1)
        {
            itkPoint = boneIt1.Value();
            pclPoint.x = itkPoint[0];
            pclPoint.y = itkPoint[1];
            pclPoint.z = itkPoint[2];
            m_data->PointsCT->points.push_back(pclPoint);
            m_data->CvPointsCT.push_back(cv::Point3d(itkPoint[0], itkPoint[1], itkPoint[2]));
        }
        ++pointsBegin;
    }

    isVTK = false;

}
*/

Registration::Registration(const vtkSmartPointer<vtkPolyData> img)
{
    m_data = new RegistrationPrivate();
    isVTK = true;

    vtkNew<vtkPolyDataConnectivityFilter> imgConnectivityFilter;
    imgConnectivityFilter->SetInputData(img);
    imgConnectivityFilter->SetExtractionModeToLargestRegion();
    imgConnectivityFilter->Update();

    vtkNew<vtkCleanPolyData> CleanImg;
    CleanImg->SetInputData(imgConnectivityFilter->GetOutput());
    CleanImg->Update();
    poly = CleanImg->GetOutput();

    /*poly = vtkSmartPointer<vtkPolyData>::New();
    poly->DeepCopy(img);*/

    vtkSmartPointer<vtkPoints> pointsList = poly->GetPoints();
    int tSize = pointsList->GetNumberOfPoints();

    for (int i = 0; i < tSize; i++)
    {
        double pnt[3];
        pointsList->GetPoint(i, pnt);
        m_data->CvPointsCT.push_back(cv::Point3d(pnt[0], pnt[1], pnt[2]));
    }
}

Registration::~Registration()
{
    m_data = NULL;
    delete m_data;
}

//RegistrationPoints::RegistrationPoints(std::vector<PointTypeITK> pPoints)
//{
//    points = pPoints;
//}

vtkSmartPointer<vtkPolyData> Registration::getContour(const vtkSmartPointer<vtkPolyData> polyData, const cv::Point3d& pNormal, const cv::Point3d& pPoint)
{
    vtkNew<vtkPlane> plane;
    vtkNew<vtkCutter> cutter;
    cutter->SetInputData(polyData);

    plane->SetNormal(pNormal.x, pNormal.y, pNormal.z);
    plane->SetOrigin(pPoint.x, pPoint.y, pPoint.z);
    cutter->SetCutFunction(plane);
    cutter->Update();

    auto contour = cutter->GetOutput();

    return contour;
}

PointTypeITK Registration::TransformMarkerPointToCtPointITK(const PointTypeITK& point) const
{
    Eigen::Vector4d tPoint;
    tPoint << point[0], point[1], point[2], 1.0;
    Eigen::Vector4d transform = m_data->mMatrix * tPoint;
    PointTypeITK result;
    result[0] = transform(0, 0);
    result[1] = transform(1, 0);
    result[2] = transform(2, 0);
    return result;
}

PointTypeITK Registration::TransformCTPointToMarkerPointITK(const PointTypeITK& point) const
{
    Eigen::Vector4d tPoint;
    tPoint << point[0], point[1], point[2], 1.0;
    Eigen::Vector4d transform = m_data->mMatrix.inverse() * tPoint;
    PointTypeITK result;
    result[0] = transform(0, 0);
    result[1] = transform(1, 0);
    result[2] = transform(2, 0);
    return result;
}


double Registration::getDistanceToBone(double x, double y, double z) const
{
    /*
    if (isVTK == false)
    {
        pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
        kdtree.setInputCloud(m_data->PointsCT);

        int K = 1;

        std::vector<int> pointIdxNKNSearch(K);
        std::vector<float> pointNKNSquaredDistance(K);
        pcl::PointXYZ refPoint(x, y, z);

        if (kdtree.nearestKSearch(refPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
        {
            return sqrt(pointNKNSquaredDistance[0]);
        }
        else
        {
            throw RegistrationExceptionCode::CAN_NOT_DETERMINE_DISTANCE_TO_3D_MESH_SURFACE;
        }
    }
    else
    {
        vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
        implicitPolyDataDistance->SetInput(poly);
        double pnt[3];
        pnt[0] = x;
        pnt[1] = y;
        pnt[2] = z;

        double myClosest[3];

        double distance = implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);
        return abs(distance);
    }
    */

    vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    implicitPolyDataDistance->SetInput(poly);
    double pnt[3];
    pnt[0] = x;
    pnt[1] = y;
    pnt[2] = z;

    double myClosest[3];

    double distance = implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);
    return abs(distance);
}

double Registration::getError() const
{
    return m_data->mError;
}

/*
std::vector<PointTypeITK> Registration::getPointsCT() const
{
    std::vector<PointTypeITK> result;
    PointTypeITK temp;
    for (int i = 0; i < m_data->PointsCT->points.size(); i++)
    {
        pcl::PointXYZ pclPoint = m_data->PointsCT->points[i];
        temp[0] = pclPoint.x;
        temp[1] = pclPoint.y;
        temp[2] = pclPoint.z;
        result.push_back(temp);
    }
    return result;
}
*/

PointTypeITK Registration::getPointInsideSphere(const PointTypeITK& center, double radius)
{
    PointTypeITK point;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distr(0.0, 1.0);
    std::uniform_real_distribution<> distr_r(0.0, radius);

    double U = distr(gen);
    double V = (distr(gen)) * 2.0 - 1.0;

    double theta = 2.0*EIGEN_PI*U;
    double phi;

    if (V <= -1.0)
    {
        phi = EIGEN_PI;
    }
    else if (V >= 1.0)
    {
        phi = 0;
    }
    else
    {
        phi = acos(V);
    }

    double R = distr_r(gen);
    point[0] = center[0] + R * sin(phi) * cos(theta);
    point[1] = center[1] + R * sin(phi) * sin(theta);
    point[2] = center[2] + R * cos(phi);
    return point;
}



PointTypeITK Registration::makeItkPoint(double x, double y, double z)
{
    PointTypeITK point;
    point[0] = x;
    point[1] = y;
    point[2] = z;
    return point;
}

itk::Rigid3DTransform<double>::Pointer Registration::getTransformMarkerToCt() const
{
    itk::Matrix< double, 3, 3 > rotation;
    itk::Vector< double, 3 > translate;

    rotation[0][0] = m_data->mMatrix(0, 0);
    rotation[0][1] = m_data->mMatrix(0, 1);
    rotation[0][2] = m_data->mMatrix(0, 2);

    rotation[1][0] = m_data->mMatrix(1, 0);
    rotation[1][1] = m_data->mMatrix(1, 1);
    rotation[1][2] = m_data->mMatrix(1, 2);

    rotation[2][0] = m_data->mMatrix(2, 0);
    rotation[2][1] = m_data->mMatrix(2, 1);
    rotation[2][2] = m_data->mMatrix(2, 2);

    translate[0] = m_data->mMatrix(0, 3);
    translate[1] = m_data->mMatrix(1, 3);
    translate[2] = m_data->mMatrix(2, 3);

    itk::Rigid3DTransform<double>::Pointer transform = itk::VersorRigid3DTransform<double>::New();
    transform->SetMatrix(rotation);
    transform->SetOffset(translate);
    return transform;
}

itk::Rigid3DTransform<double>::Pointer Registration::getTransformCtToMarker() const
{
    Eigen::Matrix4d matrix_inverse = m_data->mMatrix.inverse();
    itk::Matrix< double, 3, 3 > rotation;
    itk::Vector< double, 3 > translate;
    rotation[0][0] = matrix_inverse(0, 0);
    rotation[0][1] = matrix_inverse(0, 1);
    rotation[0][2] = matrix_inverse(0, 2);

    rotation[1][0] = matrix_inverse(1, 0);
    rotation[1][1] = matrix_inverse(1, 1);
    rotation[1][2] = matrix_inverse(1, 2);

    rotation[2][0] = matrix_inverse(2, 0);
    rotation[2][1] = matrix_inverse(2, 1);
    rotation[2][2] = matrix_inverse(2, 2);

    translate[0] = matrix_inverse(0, 3);
    translate[1] = matrix_inverse(1, 3);
    translate[2] = matrix_inverse(2, 3);

    itk::Rigid3DTransform<double>::Pointer transform = itk::VersorRigid3DTransform<double>::New();
    transform->SetMatrix(rotation);
    transform->SetOffset(translate);
    return transform;
}



void Registration::deletePointsInsideRadius(std::list<cv::Point3d>& points, const cv::Point3d& centerPoint, double radius)
{
    double distance;
    auto it = points.begin();
    while (it != points.end())
    {
        distance = RLine::getDistanceBetweenPoints(centerPoint, *it, true);
        if (distance < radius)
        {
            it = points.erase(it);
        }
        else
        {
            ++it;
        }
    }
}

int Registration::getOneVerticePosition(const std::vector<cv::Point3d>& sortPoints, const RLine& perpendicularRefL1, const RLine& perpendicularL2, double factor)
{
    RPlane planeTemp;
    planeTemp.init(perpendicularRefL1.getDirectVector(), perpendicularL2.getPoint());
    cv::Point3d interception = planeTemp.getInterceptionLinePoint(perpendicularRefL1);
    cv::Point3d V1 = perpendicularRefL1.getPoint() - interception;
    cv::Point3d V2 = perpendicularL2.getPoint() - interception;
    V1 = V1 / sqrt(V1.dot(V1));
    V2 = V2 / sqrt(V2.dot(V2));

    cv::Point3d a = interception + factor * 10000.0 * V1;
    cv::Point3d b = interception + 10000.0 * V2;

    cv::Point3d directVector = a - b;

    RLine line = RLine(directVector, a);
    double distance = -1.0;
    double temp;
    std::vector<cv::Point3d>::const_iterator it1 = sortPoints.begin();
    std::vector<cv::Point3d>::const_iterator it2 = sortPoints.end();
    int pos = -1;
    int cont = 0;
    if (sortPoints.size() < 10)
    {
        return pos;
    }

    for (; it1 != it2; ++it1)
    {
        temp = line.getDistanceFromPoint(*it1);
        if (temp > distance)
        {
            distance = temp;
            pos = cont;
        }
        cont++;
    }
    return pos;
}

//It always starts searching on the negative side until he finds the first positive evaluation.
std::vector<cv::Point3d> Registration::reduceSortPoints(const std::vector<cv::Point3d>& sortPoints, const RPlane& plane, bool isReverse)
{
    std::vector<cv::Point3d> result;

    if (isReverse == true)
    {
        for (auto it = sortPoints.rbegin(); it != sortPoints.rend(); ++it)
        {
            if (plane.eval(*it) > 0)
            {
                std::vector<cv::Point3d> newVect(it, sortPoints.rend());
                return newVect;
            }
        }
    }
    else
    {
        for (auto it = sortPoints.begin(); it != sortPoints.end(); ++it)
        {
            if (plane.eval(*it) > 0)
            {
                std::vector<cv::Point3d> newVect(it, sortPoints.end());
                return newVect;
            }
        }
    }
    return result;
}

std::vector<cv::Point3d> Registration::reduceSortPointsByAngle(const std::vector<cv::Point3d>& sortPoints, const cv::Point3d& vector, const cv::Point3d& centerPoint, double minAngle, bool isReverse)
{
    std::vector<cv::Point3d> result;
    cv::Point3d mineVector;
    double angle;

    if (isReverse == true)
    {
        for (auto it = sortPoints.rbegin(); it != sortPoints.rend(); ++it)
        {
            mineVector = *it - centerPoint;
            angle = RLine::getAngleBetweenVectors(vector, mineVector);
            if (angle < minAngle)
            {
                std::vector<cv::Point3d> newVect(it, sortPoints.rend());
                return newVect;
            }
        }
    }
    else
    {
        for (auto it = sortPoints.begin(); it != sortPoints.end(); ++it)
        {
            mineVector = *it - centerPoint;
            angle = RLine::getAngleBetweenVectors(vector, mineVector);
            if (angle < minAngle)
            {
                std::vector<cv::Point3d> newVect(it, sortPoints.end());
                return newVect;
            }
        }
    }
    return result;
}

void Registration::getPointAtRadius(const std::vector<cv::Point3d>& sortPoints, std::vector<cv::Point3d>& outPoints, const cv::Point3d& centerPoint, double squareRadius, int amount)
{
    std::list<cv::Point3d> myList(sortPoints.begin(), sortPoints.end());
    double distance, distanceTemp;
    cv::Point3d changePoint = centerPoint;
    cv::Point3d tempPoint;
    while (amount > 0)
    {
        deletePointsInsideRadius(myList, changePoint, squareRadius);
        distance = 999999;
        auto it1 = myList.begin();
        auto it2 = myList.end();
        for (; it1 != it2; ++it1)
        {
            distanceTemp = RLine::getDistanceBetweenPoints(*it1, changePoint);
            if (distanceTemp < distance)
            {
                distance = distanceTemp;
                tempPoint = *it1;
            }
        }
        outPoints.push_back(tempPoint);
        changePoint = tempPoint;
        amount--;
    }
}

cv::Mat Registration::GetTranslationRotation(const std::vector<PointTypeITK>& source, const std::vector<PointTypeITK>& target)
{
    Eigen::Matrix4d rigidTransform = m_data->InitRegistration(source, target);

    Eigen::Matrix3d rotation(3, 3);

    rotation(0, 0) = rigidTransform(0, 0);
    rotation(1, 0) = rigidTransform(1, 0);
    rotation(2, 0) = rigidTransform(2, 0);

    rotation(0, 1) = rigidTransform(0, 1);
    rotation(1, 1) = rigidTransform(1, 1);
    rotation(2, 1) = rigidTransform(2, 1);

    rotation(0, 2) = rigidTransform(0, 2);
    rotation(1, 2) = rigidTransform(1, 2);
    rotation(2, 2) = rigidTransform(2, 2);

    Eigen::Vector3d rot = rotation.eulerAngles(0, 1, 2);

    cv::Mat result(6, 1, CV_64F);
    result.at<double>(0, 0) = rigidTransform(0, 3);
    result.at<double>(1, 0) = rigidTransform(1, 3);
    result.at<double>(2, 0) = rigidTransform(2, 3);

    result.at<double>(3, 0) = rot(0);
    result.at<double>(4, 0) = rot(1);
    result.at<double>(5, 0) = rot(2);

    return result;
}

void Registration::MakeResult(const cv::Mat& pData, double pError)
{
    Eigen::AngleAxisd init_rotationZ(pData.at<double>(5, 0), Eigen::Vector3d::UnitZ());
    Eigen::AngleAxisd init_rotationY(pData.at<double>(4, 0), Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd init_rotationX(pData.at<double>(3, 0), Eigen::Vector3d::UnitX());

    Eigen::Translation3d init_translation(pData.at<double>(0, 0), pData.at<double>(1, 0), pData.at<double>(2, 0));

    Eigen::Matrix4d init_guess = (init_translation * init_rotationX * init_rotationY * init_rotationZ).matrix();

    m_data->mMatrix = init_guess;

    m_data->mError = pError;

}

itk::Rigid3DTransform<double>::Pointer Registration::GetTransformBetweenPoints(const std::vector<PointTypeITK>& source, const std::vector<PointTypeITK>& target)
{
    RegistrationPrivate myTool;

    Eigen::Matrix4d rigidTransform = myTool.InitRegistration(source, target);

    itk::Matrix< double, 3, 3 > rotation;
    itk::Vector< double, 3 > translate;
    rotation[0][0] = rigidTransform(0, 0);
    rotation[0][1] = rigidTransform(0, 1);
    rotation[0][2] = rigidTransform(0, 2);

    rotation[1][0] = rigidTransform(1, 0);
    rotation[1][1] = rigidTransform(1, 1);
    rotation[1][2] = rigidTransform(1, 2);

    rotation[2][0] = rigidTransform(2, 0);
    rotation[2][1] = rigidTransform(2, 1);
    rotation[2][2] = rigidTransform(2, 2);

    translate[0] = rigidTransform(0, 3);
    translate[1] = rigidTransform(1, 3);
    translate[2] = rigidTransform(2, 3);

    itk::Rigid3DTransform<double>::Pointer transform = itk::VersorRigid3DTransform<double>::New();
    transform->SetMatrix(rotation);
    transform->SetOffset(translate);
    return transform;
}

std::pair<cv::Point3d, double> Registration::getMinCircle(const std::vector<cv::Point3d>& pPoints, int amount, std::vector<cv::Point3d>& circlePoints)
{
    RPlane helpPlane;
    int tSize = pPoints.size();
    if (tSize == 3)
    {
        helpPlane.init(pPoints[0], pPoints[1], pPoints[2]);
    }
    else if (tSize > 3)
    {
        helpPlane = RPlane::getBestPlane(pPoints);
    }
    else
    {
        throw RegistrationExceptionCode::NOT_ENOUGH_POINTS_TO_FIT_CIRCLE;
    }

    double Z = 0;
    cv::Mat zRot = m_data->GetRotateZ(helpPlane.getNormalVector());
    std::vector<cv::Point2f> coplanar2d;

    auto it1 = pPoints.begin();
    auto it2 = pPoints.end();

    for (; it1 != it2; ++it1)
    {
        cv::Mat pointMat = zRot * m_data->cvPointToMat(*it1);
        cv::Point3d pointTemp = cv::Point3d(pointMat);
        coplanar2d.push_back(cv::Point2f(pointTemp.x, pointTemp.y));
        Z += pointTemp.z;
    }

    Z = Z / double(tSize);

    cv::Point2f center;
    float radius;

    cv::minEnclosingCircle(coplanar2d, center, radius);

    double step = (2.0 * EIGEN_PI) / double(amount);

    for (int i = 0; i < amount; i++)
    {
        double angle = double(i + 1) * step;
        double X = center.x + (radius) * cos(angle);
        double Y = center.y + (radius) * sin(angle);

        cv::Mat pointMat = zRot.inv() * m_data->cvPointToMat(cv::Point3d(X, Y, Z));
        circlePoints.push_back(cv::Point3d(pointMat));
    }

    cv::Point3d myCenter = cv::Point3d(center.x, center.y, Z);

    cv::Mat centerMat = zRot.inv() * m_data->cvPointToMat(myCenter);

    return std::make_pair(cv::Point3d(centerMat), radius);
}



/*
inline void mouseEventOccurred(const pcl::visualization::PointPickingEvent &event,
    void* viewer_void)
{
    std::cout << "Picking event active" << std::endl;
    if (event.getPointIndex() == -1)
    {
        std::cout << "Do not Point" << std::endl;
    }
    else
    {
        float x, y, z;
        event.getPoint(x, y, z);
        std::cout << x << "; " << y << "; " << z << std::endl;
    }

}

void Registration::drawCloudRGB(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud)
{

    pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
    viewer->setBackgroundColor(0, 0, 0);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(cloud);
    viewer->addPointCloud<pcl::PointXYZRGB>(cloud, rgb, "sample cloud");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "sample cloud");
    viewer->addCoordinateSystem(1.0);
    viewer->initCameraParameters();

    viewer->registerPointPickingCallback(mouseEventOccurred, (void*)viewer.get());
    while (!viewer->wasStopped())
    {
        viewer->spinOnce(100);
        std::this_thread::sleep_for(100ms);
    }
}

void Registration::drawCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::PointCloud<pcl::PointXYZ>::Ptr redPoints, bool applyCentroidToRed)
{
    std::uint8_t r = 255, g = 255, b = 255;
    std::uint32_t white = ((std::uint32_t)r << 16 | (std::uint32_t)g << 8 | (std::uint32_t)b);

    r = 255, g = 0, b = 0;
    std::uint32_t red = ((std::uint32_t)r << 16 | (std::uint32_t)g << 8 | (std::uint32_t)b);

    std::uint32_t color;

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloudRGB(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::PointXYZRGB point, centroid;
    centroid.x = 0.0;
    centroid.y = 0.0;
    centroid.z = 0.0;

    auto it11 = cloud->points.begin();
    auto it21 = cloud->points.end();
    int cont = 0;
    for (; it11 != it21; ++it11)
    {
        cont++;
        centroid.x = centroid.x + (*it11).x;
        centroid.y = centroid.y + (*it11).y;
        centroid.z = centroid.z + (*it11).z;
    }
    centroid.x = centroid.x / cont;
    centroid.y = centroid.y / cont;
    centroid.z = centroid.z / cont;

    std::cout << "Centroid: " << centroid << std::endl;

    auto it1 = cloud->points.begin();
    auto it2 = cloud->points.end();

    int redSize = redPoints->points.size();
    for (; it1 != it2; ++it1)
    {
        point.x = (*it1).x - centroid.x;
        point.y = (*it1).y - centroid.y;
        point.z = (*it1).z - centroid.z;
        color = white;


        //for (int i = 0; i < redSize; i++)
        //{
        //    //std::cout << distanceBetweenPCLPoints(point, redPoints->points[i]) << std::endl;
        //    pcl::PointXYZ redPoint = redPoints->points[i];
        //    if (applyCentroidToRed == true)
        //    {
        //        redPoint.x = redPoint.x - centroid.x;
        //        redPoint.y = redPoint.y - centroid.y;
        //        redPoint.z = redPoint.z - centroid.z;
        //    }
        //    if (distanceBetweenPCLPoints(point, redPoint) < 1.0)
        //    {
        //        color = red;
        //        break;
        //    }
        //}
        point.rgb = *reinterpret_cast<float*>(&color);
        cloudRGB->points.push_back(point);
    }

    for (int i = 0; i < redSize; i++)
    {
        pcl::PointXYZ redPoint = redPoints->points[i];
        point.x = redPoint.x - centroid.x;
        point.y = redPoint.y - centroid.y;
        point.z = redPoint.z - centroid.z;
        color = red;

        point.rgb = *reinterpret_cast<float*>(&color);
        cloudRGB->points.push_back(point);
    }

    drawCloudRGB(cloudRGB);
}

double Registration::distanceBetweenPCLPoints(pcl::PointXYZRGB& p1, pcl::PointXYZ& p2)
{
    double distance = (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z);
    return sqrt(distance);
}

///////////////////////////////////Test in header/////////////////////////////////////////////////
/*
PointTypeITK Registration::TransformPointTest(PointTypeITK& point, bool putRandom)
{
    Eigen::Vector3f newPoint;
    newPoint << point[0], point[1], point[2];
    Eigen::Vector3f tPoint = rotationTest * newPoint + translationTest;
    PointTypeITK itkPoint;
    itkPoint[0] = tPoint(0, 0);
    itkPoint[1] = tPoint(1, 0);
    itkPoint[2] = tPoint(2, 0);

    if (putRandom == true)
    {
        itkPoint[0] = itkPoint[0] + 9.0 - static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 18));
        itkPoint[1] = itkPoint[1] + 9.0 - static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 18));
        itkPoint[2] = itkPoint[2] + 9.0 - static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 18));
    }

    return itkPoint;
}

void Registration::getTestTransformPoints(std::vector<PointTypeITK>& points, std::vector<int>& pos, int amount)
{
    int size = PointsCT->points.size();
    int step = size / amount;
    int rest = size % amount;
    pcl::PointXYZ pclPoint;
    PointTypeITK itkPoint;
    for (int i = 0; i < size - rest; i += step)
    {
        pclPoint = PointsCT->points[i];
        Eigen::Vector3f point;
        point << pclPoint.x, pclPoint.y, pclPoint.z;
        Eigen::Vector3f tPoint = rotationTest * point + translationTest;

        itkPoint[0] = tPoint(0, 0);
        itkPoint[1] = tPoint(1, 0);
        itkPoint[2] = tPoint(2, 0);
        pos.push_back(i);
        points.push_back(itkPoint);
    }
}

void Registration::SetTransformMatrixTest()
{
    double x = 500.0 + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 500));
    double y = 500.0 + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 500));
    double z = 500.0 + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 500));

    double thZ = 180.0 - static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 360));
    double thY = 180.0 - static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 360));
    double thX = 180.0 - static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 360));

    Eigen::AngleAxisf init_rotationZ(thZ * M_PI / 180.0, Eigen::Vector3f::UnitZ());
    Eigen::AngleAxisf init_rotationY(thY * M_PI / 180.0, Eigen::Vector3f::UnitY());
    Eigen::AngleAxisf init_rotationX(thX * M_PI / 180.0, Eigen::Vector3f::UnitX());

    Eigen::Translation3f init_translation(x, y, z);

    Eigen::Matrix4f init_guess = (init_translation * init_rotationX * init_rotationY * init_rotationZ).matrix();

    rotationTest(0, 0) = init_guess(0, 0);
    rotationTest(0, 1) = init_guess(0, 1);
    rotationTest(0, 2) = init_guess(0, 2);
    translationTest(0, 0) = init_guess(0, 3);

    rotationTest(1, 0) = init_guess(1, 0);
    rotationTest(1, 1) = init_guess(1, 1);
    rotationTest(1, 2) = init_guess(1, 2);
    translationTest(1, 0) = init_guess(1, 3);

    rotationTest(2, 0) = init_guess(2, 0);
    rotationTest(2, 1) = init_guess(2, 1);
    rotationTest(2, 2) = init_guess(2, 2);
    translationTest(2, 0) = init_guess(2, 3);
}

void Registration::showTestResult(pcl::PointCloud<pcl::PointXYZ>::Ptr a, pcl::PointCloud<pcl::PointXYZ>::Ptr b, std::vector<int>& pos, std::string msg)
{
    std::cout << msg << "Big set size: " << a->points.size() << ", Small set size: " << b->points.size() << std::endl;
    double norm = 0;
    pcl::PointXYZ point;
    if (pos.size() == 0)
    {
        int step = a->points.size() / 20;
        for (int i = 0; i < a->points.size(); i += step)
        {
            point.x = (a->points[i].x - b->points[i].x);
            point.y = (a->points[i].y - b->points[i].y);
            point.z = (a->points[i].z - b->points[i].z);
            norm = sqrt((point.x * point.x) + (point.y * point.y) + (point.z * point.z));

            std::cout << a->points[i] << "  " << b->points[i] << "  " << norm << std::endl;
        }
    }
    else
    {
        int step = pos.size() / pos.size();
        //std::ofstream MyFile("correspondence.txt");

        for (int j = 0; j < pos.size(); j += step)
        {
            int i = pos[j];
            point.x = (a->points[i].x - b->points[j].x);
            point.y = (a->points[i].y - b->points[j].y);
            point.z = (a->points[i].z - b->points[j].z);
            norm = sqrt((point.x * point.x) + (point.y * point.y) + (point.z * point.z));

            std::cout << a->points[i] << "  " << b->points[j] << "  " << norm << std::endl;
            //MyFile << b->points[j].x << " " << b->points[j].y << " " << b->points[j].z << "  -  " << a->points[i].x << " " << a->points[i].y << " " << a->points[i].z <<std::endl;
        }
        //MyFile.close();
    }
}
*/
/////////////////////////////////////////////////////////////////////////////////

/*

void Registration::writeFile(pcl::PointCloud<pcl::PointXYZ>::Ptr points, std::string name)
{
    std::ofstream outfile(name);
    double a, b, c;
    int tSize = points->points.size();
    for (int i = 0; i < tSize; i++)
    {
        a = points->points[i].x;
        b = points->points[i].y;
        c = points->points[i].z;
        outfile << a << " " << b << " " << " " << c << "\n";
    }
}

void Registration::writeFile(std::vector<PointTypeITK> points, std::string name)
{
    std::ofstream outfile(name);
    double a, b, c;
    int tSize = points.size();
    for (int i = 0; i < tSize; i++)
    {
        a = points[i][0];
        b = points[i][1];
        c = points[i][2];
        outfile << a << " " << b << " " << " " << c << "\n";
    }
}

///////////////////////////////////////////////////////////////////

pcl::PointXYZ Registration::itkToPclPoint(PointTypeITK& itkPoint, pcl::PointXYZ& centroid)
{
    pcl::PointXYZ P1;
    P1.x = itkPoint[0] - centroid.x;
    P1.y = itkPoint[1] - centroid.y;
    P1.z = itkPoint[2] - centroid.z;
    return P1;
}

cv::Point3d Registration::pclPointToOpenCV(pcl::PointXYZ& point)
{
    cv::Point3d cvPoint;
    cvPoint.x = point.x;
    cvPoint.y = point.y;
    cvPoint.z = point.z;
    return cvPoint;
}

pcl::PointXYZ Registration::openCvToPCL(cv::Point3d& cvPoint)
{
    pcl::PointXYZ P1;
    P1.x = cvPoint.x;
    P1.y = cvPoint.y;
    P1.z = cvPoint.z;
    return P1;
}

cv::Mat Registration::pointToMat(cv::Point3d& point)
{
    cv::Mat mat(3, 1, CV_64F);
    mat.at <double>(0, 0) = point.x;
    mat.at <double>(1, 0) = point.y;
    mat.at <double>(2, 0) = point.z;
    return mat;
}

void Registration::DrawForTest(pcl::PointCloud<pcl::PointXYZ>::Ptr pointsToDrow)
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr Points;
    pcl::PointCloud<pcl::PointXYZ>::Ptr redPoints(new pcl::PointCloud<pcl::PointXYZ>);
    Points = pointsToDrow;

    cv::Point3d hipCenter(22.5, -50.5, 874.0);
    cv::Point3d kneeCenter(24.0, -64.5, 769.0);
    cv::Point3d directVector = hipCenter - kneeCenter;

    std::vector<cv::Point3d> cvPoints;
    RPlane tibiaCutPlane;
    tibiaCutPlane.init(directVector, (kneeCenter - cv::Point3d(0, 0, 5)));
    tibiaCutPlane.normalizeNormalVector();

    int tSize = Points->points.size();
    double translationZ = -(tibiaCutPlane.getBias() / tibiaCutPlane.getNormalVector().z);
    cv::Mat trasnlation(3, 1, CV_64F);
    cv::Point3d suma(0.0, 0.0, 0.0);

    for (int i = 0; i < tSize; i++)
    {
        if (tibiaCutPlane.isPointNearToPlane(pclPointToOpenCV(Points->points[i]), 1.0) == true)
        {
            cv::Point3d tempPoint = tibiaCutPlane.getProjectionPoint(pclPointToOpenCV(Points->points[i]));
            cvPoints.push_back(tempPoint);
            suma = suma + tempPoint;
        }
    }
    cv::Point3d meanP = suma / double(cvPoints.size());
    trasnlation.at <double>(0, 0) = meanP.x;
    trasnlation.at <double>(1, 0) = meanP.y;
    trasnlation.at <double>(2, 0) = meanP.z;

    ///////////////////////////////////////////////////////////////
    cv::Point3d normalXY(0.0, 0.0, 1.0);
    cv::Point3d rotationAxis = normalXY.cross(tibiaCutPlane.getNormalVector());
    rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
    double rotationAngle = RLine::getAngleBetweenVectors(normalXY, tibiaCutPlane.getNormalVector());

    cv::Mat rotation_1 = RLine::getRotateMatrix(rotationAxis, -rotationAngle);

    cv::Mat rotation_2 = RLine::getRotateMatrix(rotationAxis, rotationAngle);

    cv::Mat rotateVector_1 = rotation_1 * tibiaCutPlane.getNormalVectorMat();
    cv::Mat rotateVector_2 = rotation_2 * tibiaCutPlane.getNormalVectorMat();

    cv::Mat rotate;

    double distance_1 = RLine::getDistanceBetweenPoints(cv::Point3d(rotateVector_1), normalXY);
    double distance_2 = RLine::getDistanceBetweenPoints(cv::Point3d(rotateVector_2), normalXY);

    if (distance_1 < distance_2)
    {
        rotate = rotation_1;
    }
    else
    {
        rotate = rotation_2;
    }

    cv::Mat aa = rotate * tibiaCutPlane.getNormalVectorMat();
    double angle22 = RLine::getAngleBetweenVectors(cv::Point3d(aa), normalXY);
    std::cout << "vector: " << rotationAxis << std::endl;
    std::cout << "angle: " << rotationAngle * 180.0 / acos(-1.0) << std::endl;
    std::cout << "Rotate: " << rotate << std::endl;
    std::cout << "angle entre: : " << angle22 * 180.0 / acos(-1.0) << std::endl;


    /////////////////////////////////////////////////////////////////

    std::vector<cv::Point2f> coplanar2d;
    std::vector<cv::Point3d> coplanar3d;
    auto it11 = cvPoints.begin();
    auto it21 = cvPoints.end();
    double sumaZ = 0.0;
    ofstream MyFile("coplanar.txt");
    for (; it11 != it21; ++it11)
    {

        cv::Mat point(3, 1, CV_64F);
        point.at <double>(0, 0) = (*it11).x;
        point.at <double>(1, 0) = (*it11).y;
        point.at <double>(2, 0) = (*it11).z;
        //cv::Mat rotatePointMat = rotate * (point - trasnlation);
        cv::Mat rotatePointMat = rotate * (point);
        cv::Point3d rotatePoint(rotatePointMat.at<double>(0, 0), rotatePointMat.at<double>(1, 0), rotatePointMat.at<double>(2, 0));
        coplanar2d.push_back(cv::Point2f(float(rotatePoint.x), float(rotatePoint.y)));
        coplanar3d.push_back(rotatePoint);
        //std::cout << "Origin: " << rotatePointMat << " Change: " << rotatePoint << std::endl;
        sumaZ = sumaZ + rotatePoint.z;
        //std::cout << rotatePoint.z << std::endl;
        MyFile << rotatePoint.x << " " << rotatePoint.y << "\n";
    }

    sumaZ = sumaZ / double(cvPoints.size());
    cv::RotatedRect box = cv::fitEllipse(coplanar2d);
    cv::Point2f squareP[4];
    box.points(squareP);
    std::cout << "With: " << box.size.width << " Height:" << box.size.height << " Center: "<< box.center<< std::endl;
    std::cout << "distance: " << sqrt((squareP[0] - squareP[1]).dot(squareP[0] - squareP[1])) << std::endl;


    cv::RotatedRect rect = cv::minAreaRect(coplanar2d);
    cv::Point2f vertices[4];
    rect.points(vertices);
    std::cout << "Rectangle with: " << rect.size.width << " Height: " << rect.size.height << " Center: " << box.center << std::endl;

    for (int i = 0; i < 100; i++)
    {
        double k = double(i) / 100.0;
        cv::Point2f a = vertices[0] + k * (vertices[1] - vertices[0]);
        MyFile << a.x << " " << a.y << "\n";
        cv::Point2f b = vertices[1] + k * (vertices[2] - vertices[1]);
        MyFile << b.x << " " << b.y << "\n";
        cv::Point2f c = vertices[2] + k * (vertices[3] - vertices[2]);
        MyFile << c.x << " " << c.y << "\n";
        cv::Point2f d = vertices[3] + k * (vertices[0] - vertices[3]);
        MyFile << d.x << " " << d.y << "\n";
    }
    MyFile.close();

    cv::Point3d x_1(squareP[0].x, squareP[0].y, sumaZ);
    cv::Point3d x_2(squareP[1].x, squareP[1].y, sumaZ);
    cv::Point3d x_3(squareP[2].x, squareP[2].y, sumaZ);
    cv::Point3d x_4(squareP[3].x, squareP[3].y, sumaZ);
    cv::Point3d x_5(box.center.x, box.center.y, sumaZ);

    cv::Mat x11 = pointToMat(x_1);
    cv::Mat x22 = pointToMat(x_2);
    cv::Mat x33 = pointToMat(x_3);
    cv::Mat x44 = pointToMat(x_4);
    cv::Mat x55 = pointToMat(x_5);

    cv::Mat x1_result = rotate.inv() * x11 + trasnlation;
    cv::Mat x2_result = rotate.inv() * x22 + trasnlation;
    cv::Mat x3_result = rotate.inv() * x33 + trasnlation;
    cv::Mat x4_result = rotate.inv() * x44 + trasnlation;

    cv::Point3d x1(x1_result.at<double>(0, 0), x1_result.at<double>(1, 0), x1_result.at<double>(2, 0));
    cv::Point3d x2(x2_result.at<double>(0, 0), x2_result.at<double>(1, 0), x2_result.at<double>(2, 0));
    cv::Point3d x3(x3_result.at<double>(0, 0), x3_result.at<double>(1, 0), x3_result.at<double>(2, 0));
    cv::Point3d x4(x4_result.at<double>(0, 0), x4_result.at<double>(1, 0), x4_result.at<double>(2, 0));


    redPoints->points.push_back(openCvToPCL(x_1));
    redPoints->points.push_back(openCvToPCL(x_2));
    redPoints->points.push_back(openCvToPCL(x_3));
    redPoints->points.push_back(openCvToPCL(x_4));
    redPoints->points.push_back(openCvToPCL(x_5));

    auto it1 = coplanar3d.begin();
    auto it2 = coplanar3d.end();

    pcl::PointCloud<pcl::PointXYZ>::Ptr coplanar(new pcl::PointCloud<pcl::PointXYZ>);
    for (int i = 0; i < 100; i++)
    {
        double k = double(i) / 100.0;
        cv::Point3d a = x_1 + k * (x_2 - x_1);
        cv::Point3d b = x_2 + k * (x_3 - x_2);
        cv::Point3d c = x_3 + k * (x_4 - x_3);
        cv::Point3d d = x_4 + k * (x_1 - x_4);
        coplanar->points.push_back(openCvToPCL(a));
        coplanar->points.push_back(openCvToPCL(b));
        coplanar->points.push_back(openCvToPCL(c));
        coplanar->points.push_back(openCvToPCL(d));
    }

    int contt = 0;
    for (; it1 != it2; ++it1)
    {
        coplanar->points.push_back(openCvToPCL(*it1));
    }

    //Registration::drawCloud(coplanar, redPoints, true);

}

void Registration::generate3D(std::string path)
{
    std::vector<std::string> filesNames;
    itk::Directory::Pointer dir = itk::Directory::New();
    dir->Load(path.c_str());
    int totalDCM = dir->GetNumberOfFiles();
    std::cout << "Total IMG: " << totalDCM << std::endl;
    std::string nameIn;
    for (int i = 0; i < totalDCM; i++)
    {
        nameIn = dir->GetFile(i);
        if (nameIn.length() > 5)
        {
            nameIn = path + "\\" + nameIn;
            filesNames.push_back(nameIn);
            std::cout << nameIn << std::endl;
        }

    }

    constexpr unsigned int InputDimension = 2;
    constexpr unsigned int OutputDimension = 3;

    using PixelType = unsigned char;
    using InputImageType = itk::Image<PixelType, InputDimension>;
    using OutputImageType = itk::Image<PixelType, OutputDimension>;

    using ReaderType = itk::ImageFileReader<InputImageType>;
    ReaderType::Pointer reader = ReaderType::New();

    using FilterType = itk::TileImageFilter<InputImageType, OutputImageType>;
    FilterType::Pointer filter = FilterType::New();

    itk::FixedArray<unsigned int, OutputDimension> layout;
    layout[0] = 1;
    layout[1] = 1;
    layout[2] = 0;

    filter->SetLayout(layout);

    for (int ii = 0; ii < filesNames.size(); ++ii)
    {
        reader->SetFileName(filesNames[ii]);

        try
        {
            reader->Update();
        }
        catch (itk::ExceptionObject & e)
        {
            std::cerr << e << std::endl;
            return;
        }

        InputImageType::Pointer input = reader->GetOutput();
        input->DisconnectPipeline();

        filter->SetInput(ii, input);
    }

    constexpr PixelType defaultValue = 128;

    filter->SetDefaultPixelValue(defaultValue);

    using WriterType = itk::ImageFileWriter<OutputImageType>;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName("IMG_3D.nrrd");
    writer->SetInput(filter->GetOutput());

    try
    {
        writer->Update();
    }
    catch (itk::ExceptionObject & error)
    {
        std::cerr << "Error: " << error << std::endl;
        return;
    }
    std::cout << "Wrote the image " << std::endl << std::endl;
    return;
}*/

