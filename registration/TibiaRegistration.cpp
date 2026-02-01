#include "RLine.hpp"
#include <random>
#include "RegistrationException.hpp"
#include "TibiaRegistration.hpp"
#include "RegistrationPrivate.hpp"
#include "LeastSquaresICP.hpp"

//#include <pcl/common/common_headers.h>
//#include <pcl/registration/icp.h>

using namespace TKA::REGISTRATION;

TibiaRegistration::TibiaRegistration(const vtkSmartPointer<vtkPolyData> img, const PointTypeITK& pTibiaTubercleCT, const PointTypeITK& pLateralmalleolusCT, const PointTypeITK& pMedialmalleolusCT)
    :Registration(img)
{
    tibiaTubercleCT = pTibiaTubercleCT;
    lateralmalleolusCT = pLateralmalleolusCT;
    medialmalleolusCT = pMedialmalleolusCT;
}

TibiaRegistration::~TibiaRegistration()
{
}

bool TibiaRegistration::MakeRegistration(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pTibiaTubercleCamera, const PointTypeITK& pLateralmalleolusCamera, const PointTypeITK& pMedialmalleolusCamera, bool useRandomAlignment)
{
    std::vector<PointTypeITK> source, target;

    /*
    double near_radius = 5.0;
    PointTypeITK tubercle_temp_ct;
    vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    implicitPolyDataDistance->SetInput(poly);

    double signedDistance = m_data->GetNearPoint(implicitPolyDataDistance, tibiaTubercleCT, tubercle_temp_ct);

    if (signedDistance > near_radius)
    {
        throw RegistrationException("Check your CT landmarks. Distance from the bone model is too large.");
        return false;
    }
    */

    source.push_back(pTibiaTubercleCamera);
    source.push_back(pLateralmalleolusCamera);
    source.push_back(pMedialmalleolusCamera);

    target.push_back(tibiaTubercleCT);
    target.push_back(lateralmalleolusCT);
    target.push_back(medialmalleolusCT);

    cv::Mat data = Registration::GetTranslationRotation(source, target);

    std::vector<itk::Point<double, 3>> myBonePoints = pBonePoints;
    myBonePoints.push_back(pTibiaTubercleCamera);

    LeastSquaresICP myICP(myBonePoints);

    double error;

    if (useRandomAlignment == true)
    {
        error = myICP.LeastSquaresRandomInit(poly, data);
    }
    else
    {
		vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
		implicitPolyDataDistance->SetInput(poly);
        error = myICP.LeastSquares(implicitPolyDataDistance, data);
    }

    Registration::MakeResult(data, error);

    return true;
}

/*
bool TibiaRegistration::RegistrationLandmarksPCL(const PointTypeITK& pTibiaTubercleCamera, const PointTypeITK& pLateralmalleolusCamera, const PointTypeITK& pMedialmalleolusCamera, double& error)
{
    double near_radius = 5.0;

    std::vector<pcl::PointXYZ> knee_list_ct;
    std::vector<pcl::PointXYZ> tubercle_list_ct;
    std::vector<pcl::PointXYZ> lateral_list_ct;
    std::vector<pcl::PointXYZ> medial_list_ct;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_all;
    pcl::PointCloud<pcl::PointXYZ>::Ptr landMarkPoints(new pcl::PointCloud<pcl::PointXYZ>);

    std::vector<pcl::PointXYZ> pclPoints;

    pclPoints.push_back(m_data->itkPointToPCL(tibiaTubercleCT));
    pclPoints.push_back(m_data->itkPointToPCL(lateralmalleolusCT));
    pclPoints.push_back(m_data->itkPointToPCL(medialmalleolusCT));

    cloud_all = m_data->getNearPointsToReferencePoints(pclPoints, knee_list_ct, tubercle_list_ct, lateral_list_ct, medial_list_ct, near_radius);
    cloud_all->points.push_back(pcl::PointXYZ(lateralmalleolusCT[0], lateralmalleolusCT[1], lateralmalleolusCT[2]));
    cloud_all->points.push_back(pcl::PointXYZ(medialmalleolusCT[0], medialmalleolusCT[1], medialmalleolusCT[2]));

    std::vector<PointTypeITK> source, target;

    PointTypeITK knee_temp_ct, tubercle_temp_ct;

    if (knee_list_ct.size() > 0 && tubercle_list_ct.size() > 0)
    {
        knee_temp_ct = m_data->pclPointToItk(knee_list_ct[0]);
        tubercle_temp_ct = m_data->pclPointToItk(tubercle_list_ct[0]);
    }
    else
    {
        throw RegistrationException("Check your CT landmarks. Distance from the bone model too large.");
        return false;
    }

    source.push_back(pTibiaTubercleCamera);
    source.push_back(pLateralmalleolusCamera);
    source.push_back(pMedialmalleolusCamera);

    target.push_back(tubercle_temp_ct);
    target.push_back(lateralmalleolusCT);
    target.push_back(medialmalleolusCT);

    Eigen::Matrix4d rigidTransform = m_data->InitRegistration(source, target);

    std::vector<PointTypeITK> cameraPoints = {pTibiaTubercleCamera, pLateralmalleolusCamera, pMedialmalleolusCamera};

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

bool TibiaRegistration::RegistrationLandmarksVTK(const PointTypeITK& pTibiaTubercleCamera, const PointTypeITK& pLateralmalleolusCamera, const PointTypeITK& pMedialmalleolusCamera, double& error)
{
    vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    implicitPolyDataDistance->SetInput(poly);

    std::vector<PointTypeITK> source, target;
    PointTypeITK my_knee_ct, my_tubercle_ct, my_lateral_ct, my_medial_ct;

    //double signedDistance1 = m_data->GetNearPoint(implicitPolyDataDistance, tibiaKneeCenterCT, my_knee_ct);
    double signedDistance2 = m_data->GetNearPoint(implicitPolyDataDistance, tibiaTubercleCT, my_tubercle_ct);

    my_lateral_ct = lateralmalleolusCT;
    my_medial_ct = medialmalleolusCT;

    if (signedDistance2 > 5 )
    {
        throw RegistrationException("Check your CT landmarks. Distance from the bone model is too large.");
        return false;
    }

    source.push_back(pTibiaTubercleCamera);
    source.push_back(pLateralmalleolusCamera);
    source.push_back(pMedialmalleolusCamera);

    target.push_back(my_tubercle_ct);
    target.push_back(my_lateral_ct);
    target.push_back(my_medial_ct);

    Eigen::Matrix4d rigidTransform = m_data->InitRegistration(source, target);

    PointTypeITK temp, lateralTemp, medialTemp;

    //double a = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pTibiaKneeCenterCamera, temp);
    double b = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pTibiaTubercleCamera, temp);
    m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pLateralmalleolusCamera, lateralTemp);
    m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pMedialmalleolusCamera, medialTemp);

    double c = m_data->GetITKPointDistance(my_lateral_ct, lateralTemp);
    double d = m_data->GetITKPointDistance(my_medial_ct, medialTemp);

    std::vector<double> distances = {b, c, d};
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

bool TibiaRegistration::RegistrationLandmarks(const PointTypeITK& pTibiaKneeCenterCamera, const PointTypeITK& pTibiaTubercleCamera, const PointTypeITK& pLateralmalleolusCamera, const PointTypeITK& pMedialmalleolusCamera, double& error)
{
    if (isVTK == false)
    {
        return RegistrationLandmarksPCL(pTibiaKneeCenterCamera, pTibiaTubercleCamera, pLateralmalleolusCamera, pMedialmalleolusCamera, error);
    }
    else
    {
        return RegistrationLandmarksVTK(pTibiaKneeCenterCamera, pTibiaTubercleCamera, pLateralmalleolusCamera, pMedialmalleolusCamera, error);
    }
}

bool TibiaRegistration::MakeRegistrationLSRandomPCL(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pTibiaTubercleCamera, const PointTypeITK& pLateralmalleolusCamera, const PointTypeITK& pMedialmalleolusCamera)
{
    double tolerance = 7.0;
    double near_radius = 5.0;
    double max_distance = tolerance;

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_all(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr bone_points(new pcl::PointCloud<pcl::PointXYZ>);

    *cloud_all = *(m_data->PointsCT);
    cloud_all->points.push_back(pcl::PointXYZ(lateralmalleolusCT[0], lateralmalleolusCT[1], lateralmalleolusCT[2]));
    cloud_all->points.push_back(pcl::PointXYZ(medialmalleolusCT[0], medialmalleolusCT[1], medialmalleolusCT[2]));

    if (pBonePoints.size() == 0)
    {
        return false;
    }
    else
    {
        m_data->itkPointVectorToPCLPointVector(pBonePoints, bone_points);
    }
    bone_points->points.push_back(pcl::PointXYZ(pTibiaTubercleCamera[0], pTibiaTubercleCamera[1], pTibiaTubercleCamera[2]));

    std::vector<pcl::PointXYZ> knee_list_ct;
    std::vector<pcl::PointXYZ> tubercle_list_ct;
    std::vector<pcl::PointXYZ> lateral_list_ct;
    std::vector<pcl::PointXYZ> medial_list_ct;

    std::vector<pcl::PointXYZ> pclPoints;

    pclPoints.push_back(m_data->itkPointToPCL(tibiaTubercleCT));
    pclPoints.push_back(m_data->itkPointToPCL(lateralmalleolusCT));
    pclPoints.push_back(m_data->itkPointToPCL(medialmalleolusCT));

    m_data->getMostNearPointsToReferencePoints(pclPoints, knee_list_ct, tubercle_list_ct, lateral_list_ct, medial_list_ct, near_radius);
    Eigen::Matrix4d rigidTransform;
    std::vector<PointTypeITK> source, target;
    double radius = 2.0;

    PointTypeITK knee_temp_ct, tubercle_temp_ct;

    if (knee_list_ct.size() > 0 && tubercle_list_ct.size() > 0)
    {
        knee_temp_ct = m_data->pclPointToItk(knee_list_ct[0]);
        tubercle_temp_ct = m_data->pclPointToItk(tubercle_list_ct[0]);
    }
    else
    {
        throw RegistrationException("Check your CT landmarks. Distance from the bone model too large.");
        return false;
    }

    PointTypeITK lateral_temp_camera = pLateralmalleolusCamera;
    PointTypeITK medial_temp_camera = pMedialmalleolusCamera;

    PointTypeITK my_knee_ct = knee_temp_ct;
    PointTypeITK my_tubercle_ct = tubercle_temp_ct;
    PointTypeITK my_lateral_camera = lateral_temp_camera;
    PointTypeITK my_medial_camera = medial_temp_camera;

    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(cloud_all);

    for (int i = 0; i < 10000; i++)
    {
        source.clear();
        target.clear();

        source.push_back(pTibiaTubercleCamera);
        source.push_back(my_lateral_camera);
        source.push_back(my_medial_camera);

        target.push_back(my_tubercle_ct);
        target.push_back(lateralmalleolusCT);
        target.push_back(medialmalleolusCT);

        rigidTransform = m_data->InitRegistration(source, target);

        bone_points->points.push_back(pcl::PointXYZ(my_lateral_camera[0], my_lateral_camera[1], my_lateral_camera[2]));
        bone_points->points.push_back(pcl::PointXYZ(my_medial_camera[0], my_medial_camera[1], my_medial_camera[2]));

        pcl::PointCloud<pcl::PointXYZ>::Ptr output_cloud(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::transformPointCloud(*bone_points, *output_cloud, rigidTransform);

        bone_points->points.pop_back();
        bone_points->points.pop_back();

        float dis = m_data->getFitnessScore_max(kdtree, output_cloud);

        if (i == 0)
        {
            if (dis > tolerance)
            {
                throw RegistrationException("Check your registration points. Initial registration error too large.");
                return false;
            }
        }

        if (dis < max_distance)
        {
            m_data->mMatrix = rigidTransform;
            m_data->mError = dis;
            max_distance = dis;
            //std::cout << "Initial alignment: " << max_distance << std::endl;

            knee_temp_ct = my_knee_ct;
            tubercle_temp_ct = my_tubercle_ct;
            lateral_temp_camera = my_lateral_camera;
            medial_temp_camera = my_medial_camera;
        }

        my_knee_ct = m_data->getNearPointRandom(kdtree, cloud_all, getPointInsideSphere(knee_temp_ct, radius));
        my_tubercle_ct = m_data->getNearPointRandom(kdtree, cloud_all, getPointInsideSphere(tubercle_temp_ct, radius));
        my_lateral_camera = getPointInsideSphere(lateral_temp_camera, radius);
        my_medial_camera = getPointInsideSphere(medial_temp_camera, radius);
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
    //myBonePoints.push_back(pTibiaKneeCenterCamera);
    myBonePoints.push_back(pTibiaTubercleCamera);

    LeastSquaresICP myICP(myBonePoints);

    if (myICP.getMaxError() > m_data->mError)
    {
        myICP.setMaxError(m_data->mError);
    }

    double error = myICP.LeastSquares(m_data->PointsCT, data);

    //std::cout << "My ICP Error: " << error << std::endl;

    if (m_data->mError > error)
    {
        Registration::MakeResult(data, error);
    }

    return true;
}

bool TibiaRegistration::MakeRegistrationLSRandomVTK(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pTibiaTubercleCamera, const PointTypeITK& pLateralmalleolusCamera, const PointTypeITK& pMedialmalleolusCamera)
{
    double tolerance = 7.0;
    double near_radius = 5.0;
    double max_distance = tolerance;

    Eigen::Matrix4d rigidTransform;
    std::vector<PointTypeITK> source, target;
    double radius = 2.0;

    PointTypeITK knee_temp_ct, tubercle_temp_ct;

    vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    implicitPolyDataDistance->SetInput(Registration::poly);

    //double signedDistance1 = m_data->GetNearPoint(implicitPolyDataDistance, tibiaKneeCenterCT, knee_temp_ct);
    double signedDistance2 = m_data->GetNearPoint(implicitPolyDataDistance, tibiaTubercleCT, tubercle_temp_ct);

    if (signedDistance2 > 5)
    {
        throw RegistrationException("Check your CT landmarks. Distance from the bone model is too large.");
        return false;
    }

    PointTypeITK lateral_temp_camera = pLateralmalleolusCamera;
    PointTypeITK medial_temp_camera = pMedialmalleolusCamera;

    PointTypeITK my_knee_ct = knee_temp_ct;
    PointTypeITK my_tubercle_ct = tubercle_temp_ct;
    PointTypeITK my_lateral_camera = lateral_temp_camera;
    PointTypeITK my_medial_camera = medial_temp_camera;

    std::vector<Eigen::Vector4d> bone_points;
    m_data->itkPointVectorToEigenPointVector4(pBonePoints, bone_points);

    for (int i = 0; i < 10000; i++)
    {
        source.clear();
        target.clear();

        //source.push_back(pTibiaKneeCenterCamera);
        source.push_back(pTibiaTubercleCamera);
        source.push_back(my_lateral_camera);
        source.push_back(my_medial_camera);

        //target.push_back(my_knee_ct);
        target.push_back(my_tubercle_ct);
        target.push_back(lateralmalleolusCT);
        target.push_back(medialmalleolusCT);

        rigidTransform = m_data->InitRegistration(source, target);

        Eigen::Vector4d pointMed, pointLat;
        pointLat << my_lateral_camera[0], my_lateral_camera[1], my_lateral_camera[2], 1;
        pointMed << my_medial_camera[0], my_medial_camera[1], my_medial_camera[2], 1;
        bone_points.push_back(pointLat);
        bone_points.push_back(pointMed);

        std::vector<Eigen::Vector4d> output_cloud;
        m_data->transformVectorEigen4(bone_points, output_cloud, rigidTransform);

        Eigen::Vector4d transformMed = output_cloud[output_cloud.size() - 1];
        Eigen::Vector4d transformLat = output_cloud[output_cloud.size() - 2];

        Eigen::Vector4d medCT(medialmalleolusCT[0], medialmalleolusCT[1], medialmalleolusCT[2], 1);
        Eigen::Vector4d latCT(lateralmalleolusCT[0], lateralmalleolusCT[1], lateralmalleolusCT[2], 1);

        Eigen::Vector4d diffMed = transformMed - medCT;
        Eigen::Vector4d diffLat = transformLat - latCT;

        bone_points.pop_back();
        bone_points.pop_back();
        output_cloud.pop_back();
        output_cloud.pop_back();

        float disMed = sqrt(diffMed.dot(diffMed));
        float disLat = sqrt(diffLat.dot(diffLat));
        float dis = m_data->getFitnessScore_max(implicitPolyDataDistance, output_cloud);

        dis = std::max({ dis, disMed, disLat });

        if (i == 0)
        {
            if (dis > tolerance)
            {
                throw RegistrationException("Check your registration points. Initial registration error too large.");
                return false;
            }
        }

        if (dis < max_distance)
        {
            m_data->mMatrix = rigidTransform;
            m_data->mError = dis;
            max_distance = dis;
            //std::cout << "Initial alignment: " << max_distance << std::endl;

            knee_temp_ct = my_knee_ct;
            tubercle_temp_ct = my_tubercle_ct;
            lateral_temp_camera = my_lateral_camera;
            medial_temp_camera = my_medial_camera;
        }

        m_data->GetNearPoint(implicitPolyDataDistance, getPointInsideSphere(knee_temp_ct, radius), my_knee_ct);
        m_data->GetNearPoint(implicitPolyDataDistance, getPointInsideSphere(tubercle_temp_ct, radius), my_tubercle_ct);
        my_lateral_camera = getPointInsideSphere(lateral_temp_camera, radius);
        my_medial_camera = getPointInsideSphere(medial_temp_camera, radius);
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
    //myBonePoints.push_back(pTibiaKneeCenterCamera);
    myBonePoints.push_back(pTibiaTubercleCamera);

    LeastSquaresICP myICP(myBonePoints);

    if (myICP.getMaxError() > m_data->mError)
    {
        myICP.setMaxError(m_data->mError);
    }

    double error = myICP.LeastSquares(poly, data);

    if (m_data->mError > error)
    {
        Registration::MakeResult(data, error);
    }

    return true;
}
*/