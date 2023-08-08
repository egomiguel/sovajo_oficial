#include "RLine.hpp"
#include "RPlane.hpp"
#include <fstream>
//#include <pcl/filters/statistical_outlier_removal.h>
//#include <pcl/recognition/ransac_based/trimmed_icp.h>
//#include <random>
//#include <chrono>
#include "RegistrationException.hpp"
#include "FemurRegistration.hpp"
#include "RegistrationPrivate.hpp"
//#include <pcl/registration/icp.h>
//#include <pcl/exceptions.h>
#include "LeastSquaresICP.hpp"
#include "LeastSquaresPoly.hpp"
//#include "CoherentPoint/CoherentPointDrift.hpp"

using namespace PKA::REGISTRATION;

FemurRegistration::FemurRegistration(const vtkSmartPointer<vtkPolyData> img, const PointTypeITK& pHipCenterCT, const PointTypeITK& pKneeCenterCT, const PointTypeITK& pEpicondyleCT)
    :Registration(img)
{
    hipCenterCT = pHipCenterCT;
    kneeCenterCT = pKneeCenterCT;
    epicondyleCT = pEpicondyleCT;
}

FemurRegistration::~FemurRegistration()
{

}

bool FemurRegistration::MakeRegistration(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pHipCamera, const PointTypeITK& pKneeCenterCamera, const PointTypeITK& pEpicondyleCamera, bool useRandomAlignment)
{
    std::vector<PointTypeITK> source, target;

    /*
    double near_radius = 5.0;
    PointTypeITK my_knee_ct, my_medial_ct;
    vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    implicitPolyDataDistance->SetInput(Registration::poly);

    double signedDistance1 = m_data->GetNearPoint(implicitPolyDataDistance, kneeCenterCT, my_knee_ct);
    double signedDistance2 = m_data->GetNearPoint(implicitPolyDataDistance, medialEpicondyleCT, my_medial_ct);

    if (signedDistance1 > near_radius || signedDistance2 > near_radius)
    {
        throw RegistrationException("Check your CT landmarks. Distance from the bone model is too large.");
        return false;
    }*/

    source.push_back(pKneeCenterCamera);
    source.push_back(pEpicondyleCamera);
    source.push_back(pHipCamera);

    target.push_back(kneeCenterCT);
    target.push_back(epicondyleCT);
    target.push_back(hipCenterCT);

    cv::Mat data = Registration::GetTranslationRotation(source, target);

    std::vector<itk::Point<double, 3>> myBonePoints = pBonePoints;
    myBonePoints.push_back(pKneeCenterCamera);
    myBonePoints.push_back(pEpicondyleCamera);

    LeastSquaresICP myICP(myBonePoints);

    double error;

    if (useRandomAlignment == true)
    {
        error = myICP.LeastSquaresRandomInit(Registration::poly, data);
    }
    else
    {
        error = myICP.LeastSquares(Registration::poly, data);
    }

    Registration::MakeResult(data, error);

    return true;
}


/*
bool FemurRegistration::RegistrationLandmarksPCL(const PointTypeITK& pHipCamera, const PointTypeITK& pCortexCamera, const PointTypeITK& pKneeCenterCamera, const PointTypeITK& pLateralEpicondyleCamera, const PointTypeITK& pMedialEpicondyleCamera, double& error)
{
    double near_radius = 5.0;

    std::vector<PointTypeITK> source, target;
    PointTypeITK cortex_temp_ct, knee_temp_ct, lateral_temp_ct, medial_temp_ct;

    std::vector<pcl::PointXYZ> cortex_list_ct;
    std::vector<pcl::PointXYZ> knee_list_ct;
    std::vector<pcl::PointXYZ> lateralEpi_list_ct;
    std::vector<pcl::PointXYZ> medialEpi_list_ct;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_all;
    pcl::PointCloud<pcl::PointXYZ>::Ptr landMarkPoints(new pcl::PointCloud<pcl::PointXYZ>);

    std::vector<pcl::PointXYZ> pclPoints;

    pclPoints.push_back(m_data->itkPointToPCL(cortexCT));
    pclPoints.push_back(m_data->itkPointToPCL(kneeCenterCT));
    pclPoints.push_back(m_data->itkPointToPCL(lateralEpicondyleCT));
    pclPoints.push_back(m_data->itkPointToPCL(medialEpicondyleCT));

    cloud_all = m_data->getNearPointsToReferencePoints(pclPoints, cortex_list_ct, knee_list_ct, lateralEpi_list_ct, medialEpi_list_ct, near_radius);
    cloud_all->points.push_back(pcl::PointXYZ(hipCenterCT[0], hipCenterCT[1], hipCenterCT[2]));

    if (cortex_list_ct.size() > 0 && knee_list_ct.size() > 0 && lateralEpi_list_ct.size() > 0 && medialEpi_list_ct.size() > 0)
    {
        cortex_temp_ct = m_data->pclPointToItk(cortex_list_ct[0]);
        knee_temp_ct = m_data->pclPointToItk(knee_list_ct[0]);
        lateral_temp_ct = m_data->pclPointToItk(lateralEpi_list_ct[0]);
        medial_temp_ct = m_data->pclPointToItk(medialEpi_list_ct[0]);
    }
    else
    {
        throw RegistrationException("Check your CT landmarks. Distance from the bone model is too large.");
        return false;
    }

    source.push_back(pCortexCamera);
    source.push_back(pKneeCenterCamera);
    source.push_back(pLateralEpicondyleCamera);
    source.push_back(pMedialEpicondyleCamera);
    source.push_back(pHipCamera);

    target.push_back(cortex_temp_ct);
    target.push_back(knee_temp_ct);
    target.push_back(lateral_temp_ct);
    target.push_back(medial_temp_ct);
    target.push_back(hipCenterCT);

    Eigen::Matrix4d rigidTransform = m_data->InitRegistration(source, target);

    std::vector<PointTypeITK> cameraPoints = { pCortexCamera, pKneeCenterCamera, pLateralEpicondyleCamera, pMedialEpicondyleCamera, pHipCamera };

    m_data->itkPointVectorToPCLPointVector(cameraPoints, landMarkPoints);

    pcl::PointCloud<pcl::PointXYZ>::Ptr output_cloud(new pcl::PointCloud<pcl::PointXYZ>);

    pcl::transformPointCloud(*landMarkPoints, *output_cloud, rigidTransform);

    error = m_data->getFitnessScore_max(cloud_all, output_cloud);

    if (error > 7)
    {
        return false;
    }
    else
    {
        return true;
    }
}


bool FemurRegistration::RegistrationLandmarksVTK(const PointTypeITK& pHipCamera, const PointTypeITK& pCortexCamera, const PointTypeITK& pKneeCenterCamera, const PointTypeITK& pLateralEpicondyleCamera, const PointTypeITK& pMedialEpicondyleCamera, double& error)
{
    vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    implicitPolyDataDistance->SetInput(Registration::poly);

    std::vector<PointTypeITK> source, target;
    PointTypeITK my_cortex_ct, my_knee_ct, my_lateral_ct, my_medial_ct;

    double signedDistance1 = m_data->GetNearPoint(implicitPolyDataDistance, cortexCT, my_cortex_ct);
    double signedDistance2 = m_data->GetNearPoint(implicitPolyDataDistance, kneeCenterCT, my_knee_ct);
    double signedDistance3 = m_data->GetNearPoint(implicitPolyDataDistance, lateralEpicondyleCT, my_lateral_ct);
    double signedDistance4 = m_data->GetNearPoint(implicitPolyDataDistance, medialEpicondyleCT, my_medial_ct);

    if (signedDistance1 > 5 || signedDistance2 > 5 || signedDistance3 > 5 || signedDistance4 > 5)
    {
        throw RegistrationException("Check your CT landmarks. Distance from the bone model is too large.");
        return false;
    }

    source.push_back(pCortexCamera);
    source.push_back(pKneeCenterCamera);
    source.push_back(pLateralEpicondyleCamera);
    source.push_back(pMedialEpicondyleCamera);
    source.push_back(pHipCamera);

    target.push_back(my_cortex_ct);
    target.push_back(my_knee_ct);
    target.push_back(my_lateral_ct);
    target.push_back(my_medial_ct);
    target.push_back(hipCenterCT);

    Eigen::Matrix4d rigidTransform = m_data->InitRegistration(source, target);

    PointTypeITK temp;

    double a = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pCortexCamera, temp);
    double b = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pKneeCenterCamera, temp);
    double c = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pLateralEpicondyleCamera, temp);
    double d = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pMedialEpicondyleCamera, temp);
    m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pHipCamera, temp);
    double e = m_data->GetITKPointDistance(hipCenterCT, temp);

    std::vector<double> distances = { a, b, c, d, e };
    std::sort(distances.rbegin(), distances.rend());

    error = distances[0];

    if (error > 7)
    {
        return false;
    }
    else
    {
        return true;
    }
}


bool FemurRegistration::RegistrationLandmarks(const PointTypeITK& pHipCamera, const PointTypeITK& pCortexCamera, const PointTypeITK& pKneeCenterCamera, const PointTypeITK& pLateralEpicondyleCamera, const PointTypeITK& pMedialEpicondyleCamera, double& error)
{
    if (isVTK == true)
    {
        return RegistrationLandmarksVTK(pHipCamera, pCortexCamera, pKneeCenterCamera, pLateralEpicondyleCamera, pMedialEpicondyleCamera, error);
    }
    else
    {
        return RegistrationLandmarksPCL(pHipCamera, pCortexCamera, pKneeCenterCamera, pLateralEpicondyleCamera, pMedialEpicondyleCamera, error);
    }
}

bool FemurRegistration::MakeRegistrationLSRandomPCL(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pHipCamera, const PointTypeITK& pKneeCenterCamera, const PointTypeITK& pMedialEpicondyleCamera)
{
    double tolerance = 7.0;
    double near_radius = 5.0;
    double max_distance = tolerance;

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_all(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr bone_points(new pcl::PointCloud<pcl::PointXYZ>);

    (*cloud_all) = *(m_data->PointsCT);
    cloud_all->points.push_back(pcl::PointXYZ(hipCenterCT[0], hipCenterCT[1], hipCenterCT[2]));

    if (pBonePoints.size() == 0)
    {
        return false;
    }
    else
    {
        m_data->itkPointVectorToPCLPointVector(pBonePoints, bone_points);
    }

    //bone_points->points.push_back(pcl::PointXYZ(pLateralEpicondyleCamera[0], pLateralEpicondyleCamera[1], pLateralEpicondyleCamera[2]));
    bone_points->points.push_back(pcl::PointXYZ(pMedialEpicondyleCamera[0], pMedialEpicondyleCamera[1], pMedialEpicondyleCamera[2]));

    std::vector<pcl::PointXYZ> knee_list_ct;
    std::vector<pcl::PointXYZ> lateralEpi_list_ct;
    std::vector<pcl::PointXYZ> medialEpi_list_ct;

    std::vector<pcl::PointXYZ> pclPoints;

    //pclPoints.push_back(m_data->itkPointToPCL(cortexCT));
    pclPoints.push_back(m_data->itkPointToPCL(kneeCenterCT));
    //pclPoints.push_back(m_data->itkPointToPCL(lateralEpicondyleCT));
    pclPoints.push_back(m_data->itkPointToPCL(medialEpicondyleCT));

    m_data->getMostNearPointsToReferencePoints(pclPoints, cortex_list_ct, knee_list_ct, lateralEpi_list_ct, medialEpi_list_ct, near_radius);

    Eigen::Matrix4d rigidTransform;

    std::vector<PointTypeITK> source, target;

    double radius = 2.0;

    PointTypeITK cortex_temp_ct, knee_temp_ct, lateral_temp_ct, medial_temp_ct;

    if (cortex_list_ct.size() > 0 && knee_list_ct.size() > 0 && lateralEpi_list_ct.size() > 0 && medialEpi_list_ct.size() > 0)
    {
        cortex_temp_ct = m_data->pclPointToItk(cortex_list_ct[0]);
        knee_temp_ct = m_data->pclPointToItk(knee_list_ct[0]);
        lateral_temp_ct = m_data->pclPointToItk(lateralEpi_list_ct[0]);
        medial_temp_ct = m_data->pclPointToItk(medialEpi_list_ct[0]);
    }
    else
    {
        throw RegistrationException("Check your CT landmarks. Distance from the bone model is too large.");
        return false;
    }

    PointTypeITK my_cortex_ct = cortex_temp_ct;
    PointTypeITK my_knee_ct = knee_temp_ct;
    PointTypeITK my_lateral_ct = lateral_temp_ct;
    PointTypeITK my_medial_ct = medial_temp_ct;

    PointTypeITK hip_temp = pHipCamera;
    PointTypeITK hip_camera = pHipCamera;

    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(cloud_all);

    for (int i = 0; i < 10000; i++)
    {
        source.clear();
        target.clear();

        //source.push_back(pCortexCamera);
        source.push_back(pKneeCenterCamera);
        //source.push_back(pLateralEpicondyleCamera);
        source.push_back(pMedialEpicondyleCamera);
        source.push_back(hip_temp);

        //target.push_back(my_cortex_ct);
        target.push_back(my_knee_ct);
        //target.push_back(my_lateral_ct);
        target.push_back(my_medial_ct);
        target.push_back(hipCenterCT);

        bone_points->points.push_back(pcl::PointXYZ(hip_temp[0], hip_temp[1], hip_temp[2]));

        rigidTransform = m_data->InitRegistration(source, target);

        pcl::PointCloud<pcl::PointXYZ>::Ptr output_cloud(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::transformPointCloud(*bone_points, *output_cloud, rigidTransform);

        bone_points->points.pop_back();

        float dis = m_data->getFitnessScore_max(kdtree, output_cloud);

        if (i == 0)
        {
            if (dis > tolerance)
            {
                throw RegistrationException("Check your registration points. Initial registration error is too large.");
                return false;
            }
        }

        if (dis < max_distance)
        {
            m_data->mMatrix = rigidTransform;
            m_data->mError = dis;
            max_distance = dis;
            //std::cout << "Initial alignment: " << max_distance << std::endl;
            cortex_temp_ct = my_cortex_ct;
            knee_temp_ct = my_knee_ct;
            lateral_temp_ct = my_lateral_ct;
            medial_temp_ct = my_medial_ct;

            hip_camera = hip_temp;
        }

        hip_temp = getPointInsideSphere(hip_camera, radius);
        my_cortex_ct = m_data->getNearPointRandom(kdtree, cloud_all, getPointInsideSphere(cortex_temp_ct, radius));
        my_knee_ct = m_data->getNearPointRandom(kdtree, cloud_all, getPointInsideSphere(knee_temp_ct, radius));
        my_lateral_ct = m_data->getNearPointRandom(kdtree, cloud_all, getPointInsideSphere(lateral_temp_ct, radius));
        my_medial_ct = m_data->getNearPointRandom(kdtree, cloud_all, getPointInsideSphere(medial_temp_ct, radius));
    }

    Eigen::Matrix3d rotation(3, 3);

    rotation(0, 0) = m_data->mMatrix(0, 0);
    rotation(1, 0) = m_data->mMatrix(1, 0);
    rotation(2, 0) = m_data->mMatrix(2, 0);

    rotation(0, 1) = m_data->mMatrix(0, 1);
    rotation(1, 1) = m_data->mMatrix(1, 1);
    rotation(2, 1) = m_data->mMatrix(2, 1);

    rotation(0, 2) = m_data->mMatrix(0, 2);
    rotation(1, 2) = m_data->mMatrix(1, 2);
    rotation(2, 2) = m_data->mMatrix(2, 2);

    Eigen::Vector3d angles = rotation.eulerAngles(0, 1, 2);

    cv::Mat data(6, 1, CV_64F);
    data.at<double>(0, 0) = m_data->mMatrix(0, 3);
    data.at<double>(1, 0) = m_data->mMatrix(1, 3);
    data.at<double>(2, 0) = m_data->mMatrix(2, 3);

    data.at<double>(3, 0) = angles(0);
    data.at<double>(4, 0) = angles(1);
    data.at<double>(5, 0) = angles(2);

    std::vector<itk::Point<double, 3>> myBonePoints = pBonePoints;
    //myBonePoints.push_back(pCortexCamera);
    //myBonePoints.push_back(pKneeCenterCamera);
    //myBonePoints.push_back(pLateralEpicondyleCamera);
    //myBonePoints.push_back(pMedialEpicondyleCamera);

    LeastSquaresICP myICP(myBonePoints);

    if (myICP.getMaxError() > m_data->mError)
    {
        myICP.setMaxError(m_data->mError);
    }

    double error = myICP.LeastSquares(m_data->PointsCT, data);

    if (m_data->mError > error)
    {
        Registration::MakeResult(data, error);
    }

    return true;
}

bool FemurRegistration::MakeRegistrationLSRandomVTK(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pHipCamera, const PointTypeITK& pKneeCenterCamera, const PointTypeITK& pMedialEpicondyleCamera)
{
    double tolerance = 7.0;
    double near_radius = 5.0;
    double max_distance = tolerance;

    Eigen::Matrix4d rigidTransform;

    std::vector<PointTypeITK> source, target;

    double radius = 2.0;

    PointTypeITK my_cortex_ct, my_knee_ct, my_lateral_ct, my_medial_ct;
    PointTypeITK cortex_temp_ct, knee_temp_ct, lateral_temp_ct, medial_temp_ct;

    vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    implicitPolyDataDistance->SetInput(Registration::poly);

    //double signedDistance1 = m_data->GetNearPoint(implicitPolyDataDistance, cortexCT, my_cortex_ct);
    double signedDistance2 = m_data->GetNearPoint(implicitPolyDataDistance, kneeCenterCT, my_knee_ct);
    //double signedDistance3 = m_data->GetNearPoint(implicitPolyDataDistance, lateralEpicondyleCT, my_lateral_ct);
    double signedDistance4 = m_data->GetNearPoint(implicitPolyDataDistance, medialEpicondyleCT, my_medial_ct);

    if (signedDistance2 > 5 || signedDistance4 > 5)
    {
        throw RegistrationException("Check your CT landmarks. Distance from the bone model is too large.");
        return false;
    }

    PointTypeITK hip_temp = pHipCamera;
    PointTypeITK hip_camera = pHipCamera;

    std::vector<Eigen::Vector4d> bone_points;

    m_data->itkPointVectorToEigenPointVector4(pBonePoints, bone_points);

    for (int i = 0; i < 10000; i++)
    {
        source.clear();
        target.clear();

        //source.push_back(pCortexCamera);
        source.push_back(pKneeCenterCamera);
        //source.push_back(pLateralEpicondyleCamera);
        source.push_back(pMedialEpicondyleCamera);
        source.push_back(hip_temp);

        //target.push_back(my_cortex_ct);
        target.push_back(my_knee_ct);
        //target.push_back(my_lateral_ct);
        target.push_back(my_medial_ct);
        target.push_back(hipCenterCT);

        Eigen::Vector4d newPoint;
        newPoint << hip_temp[0], hip_temp[1], hip_temp[2], 1;
        bone_points.push_back(newPoint);

        rigidTransform = m_data->InitRegistration(source, target);

        std::vector<Eigen::Vector4d> output_cloud;
        m_data->transformVectorEigen4(bone_points, output_cloud, rigidTransform);

        Eigen::Vector4d transformHipTemp = output_cloud[output_cloud.size() - 1];
        Eigen::Vector4d hipCT(hipCenterCT[0], hipCenterCT[1], hipCenterCT[2], 1);
        Eigen::Vector4d diff = transformHipTemp - hipCT;

        output_cloud.pop_back();
        bone_points.pop_back();

        float disHip = sqrt(diff.dot(diff));
        float dis = m_data->getFitnessScore_max(implicitPolyDataDistance, output_cloud);

        if (dis < disHip)
        {
            dis = disHip;
        }

        if (i == 0)
        {
            if (dis > tolerance)
            {
                throw RegistrationException("Check your registration points. Initial registration error is too large.");
                return false;
            }
        }

        if (dis < max_distance)
        {
            m_data->mMatrix = rigidTransform;
            m_data->mError = dis;
            max_distance = dis;
            //std::cout << "Initial alignment: " << max_distance << std::endl;
            cortex_temp_ct = my_cortex_ct;
            knee_temp_ct = my_knee_ct;
            lateral_temp_ct = my_lateral_ct;
            medial_temp_ct = my_medial_ct;

            hip_camera = hip_temp;
        }

        hip_temp = getPointInsideSphere(hip_camera, radius);

        m_data->GetNearPoint(implicitPolyDataDistance, getPointInsideSphere(cortex_temp_ct, radius), my_cortex_ct);
        m_data->GetNearPoint(implicitPolyDataDistance, getPointInsideSphere(knee_temp_ct, radius), my_knee_ct);
        m_data->GetNearPoint(implicitPolyDataDistance, getPointInsideSphere(lateral_temp_ct, radius), my_lateral_ct);
        m_data->GetNearPoint(implicitPolyDataDistance, getPointInsideSphere(medial_temp_ct, radius), my_medial_ct);
    }

    Eigen::Matrix3d rotation(3, 3);

    rotation(0, 0) = m_data->mMatrix(0, 0);
    rotation(1, 0) = m_data->mMatrix(1, 0);
    rotation(2, 0) = m_data->mMatrix(2, 0);

    rotation(0, 1) = m_data->mMatrix(0, 1);
    rotation(1, 1) = m_data->mMatrix(1, 1);
    rotation(2, 1) = m_data->mMatrix(2, 1);

    rotation(0, 2) = m_data->mMatrix(0, 2);
    rotation(1, 2) = m_data->mMatrix(1, 2);
    rotation(2, 2) = m_data->mMatrix(2, 2);

    Eigen::Vector3d angles = rotation.eulerAngles(0, 1, 2);

    cv::Mat data(6, 1, CV_64F);
    data.at<double>(0, 0) = m_data->mMatrix(0, 3);
    data.at<double>(1, 0) = m_data->mMatrix(1, 3);
    data.at<double>(2, 0) = m_data->mMatrix(2, 3);

    data.at<double>(3, 0) = angles(0);
    data.at<double>(4, 0) = angles(1);
    data.at<double>(5, 0) = angles(2);

    std::vector<itk::Point<double, 3>> myBonePoints = pBonePoints;
    //myBonePoints.push_back(pCortexCamera);
    myBonePoints.push_back(pKneeCenterCamera);
    //myBonePoints.push_back(pLateralEpicondyleCamera);
    myBonePoints.push_back(pMedialEpicondyleCamera);

    LeastSquaresICP myICP(myBonePoints);

    if (myICP.getMaxError() > m_data->mError)
    {
        myICP.setMaxError(m_data->mError);
    }

    double error = myICP.LeastSquares(Registration::poly, data);

    if (m_data->mError > error)
    {
        Registration::MakeResult(data, error);
    }

    return true;
}

*/