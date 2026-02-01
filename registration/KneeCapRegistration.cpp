#include "KneeCapRegistration.hpp"
#include "RegistrationException.hpp"
#include "RegistrationPrivate.hpp"
//#include <pcl/registration/icp.h>
//#include <pcl/exceptions.h>
#include "LeastSquaresICP.hpp"
#include "opencv2/imgproc.hpp"
//#include "CoherentPoint/CoherentPointDrift.hpp"
#include "vtkKdTreePointLocator.h"

using namespace TKA::REGISTRATION;

/*
KneeCapRegistration::KneeCapRegistration(const std::vector<RegistrationImageType::Pointer>& img, const PointTypeITK& pHipCenterCT, const PointTypeITK& pKneeCenterCT, const PointTypeITK& pLateralEpiCT, const PointTypeITK& pMedialEpiCT)
    :Registration(img)
{
    fitPlane = RPlane::getBestPlane(m_data->CvPointsCT);

    if (fitPlane.getIsInit() == false)
    {
        throw RegistrationExceptionCode::CAN_NOT_DEFINE_FIT_PLANE_TO_PATELLA;
    }

    cv::Point3d axis = m_data->itkPointToCV(pHipCenterCT) - m_data->itkPointToCV(pKneeCenterCT);

    RPlane planeHelp;
    planeHelp.init(axis, m_data->itkPointToCV(pLateralEpiCT));

    cv::Point3d  MedialEpiCT_Perp = planeHelp.getProjectionPoint(m_data->itkPointToCV(pMedialEpiCT));

    cv::Point3d TEA = m_data->itkPointToCV(pLateralEpiCT) - MedialEpiCT_Perp;

    medLatVector = fitPlane.getProjectionVector(TEA);

    //RPlane sagital;
    //sagital.init(medLatVector, m_data->itkPointToCV(pKneeCenterCT));

    //cv::Point3d hipProj = sagital.getProjectionPoint(m_data->itkPointToCV(pHipCenterCT));
    //cv::Point3d newAxis = hipProj - m_data->itkPointToCV(pKneeCenterCT);
    cv::Point3d newAxis = m_data->itkPointToCV(pHipCenterCT) - m_data->itkPointToCV(pKneeCenterCT);

    kneeHipVector = fitPlane.getProjectionVector(newAxis);

    kneeHipVector = kneeHipVector / sqrt(kneeHipVector.dot(kneeHipVector));
    medLatVector = medLatVector / sqrt(medLatVector.dot(medLatVector));

    transformZ = GetRotateZ(fitPlane.getNormalVector());

    GetMainPoints();
}
*/

KneeCapRegistration::KneeCapRegistration(const vtkSmartPointer<vtkPolyData> img, const PointTypeITK& pHipCenterCT, const PointTypeITK& pKneeCenterCT, const PointTypeITK& pLateralEpiCT, const PointTypeITK& pMedialEpiCT)
    : Registration(img)
{
    fitPlane = RPlane::getBestPlane(m_data->CvPointsCT);

    if (fitPlane.getIsInit() == false)
    {
        throw RegistrationExceptionCode::CAN_NOT_DEFINE_FIT_PLANE_TO_PATELLA;
    }

    cv::Point3d axis = m_data->itkPointToCV(pHipCenterCT) - m_data->itkPointToCV(pKneeCenterCT);

    RPlane planeHelp;
    planeHelp.init(axis, m_data->itkPointToCV(pLateralEpiCT));

    cv::Point3d  MedialEpiCT_Perp = planeHelp.getProjectionPoint(m_data->itkPointToCV(pMedialEpiCT));

    cv::Point3d TEA = m_data->itkPointToCV(pLateralEpiCT) - MedialEpiCT_Perp;

    medLatVector = fitPlane.getProjectionVector(TEA);

    /*RPlane sagital;
    sagital.init(medLatVector, m_data->itkPointToCV(pKneeCenterCT));

    cv::Point3d hipProj = sagital.getProjectionPoint(m_data->itkPointToCV(pHipCenterCT));*/
    cv::Point3d newAxis = m_data->itkPointToCV(pHipCenterCT) - m_data->itkPointToCV(pKneeCenterCT);

    kneeHipVector = fitPlane.getProjectionVector(newAxis);

    kneeHipVector = kneeHipVector / sqrt(kneeHipVector.dot(kneeHipVector));
    medLatVector = medLatVector / sqrt(medLatVector.dot(medLatVector));

    transformZ = GetRotateZ(fitPlane.getNormalVector());

    GetMainPoints();
}

KneeCapRegistration::~KneeCapRegistration()
{
}

void KneeCapRegistration::GetMainPoints()
{
    std::vector<cv::Point3d> circlePoints;
    std::vector<cv::Point2f> transformPoints;
    std::vector<RLine> myLines;
    vtkNew<vtkPoints> projPoints;

    double z = 0;

    std::vector<cv::Point3d>::iterator it1, it2;
    it1 = m_data->CvPointsCT.begin();
    it2 = m_data->CvPointsCT.end();

    for (; it1 != it2; ++it1)
    {
        cv::Point3d temp = fitPlane.getProjectionPoint(*it1);
        double pnt[3];
        pnt[0] = temp.x;
        pnt[1] = temp.y;
        pnt[2] = temp.z;

        projPoints->InsertNextPoint(pnt);

        cv::Mat tempMat(3, 1, CV_64F);
        tempMat.at<double>(0, 0) = temp.x;
        tempMat.at<double>(1, 0) = temp.y;
        tempMat.at<double>(2, 0) = temp.z;

        cv::Mat tempTransMat = transformZ * tempMat;
        cv::Point3d tempTrans = cv::Point3d(tempTransMat);
        z += tempTrans.z;
        transformPoints.push_back(cv::Point2f(tempTrans.x, tempTrans.y));
    }

    if (m_data->CvPointsCT.size() > 0)
    {
        z = z / double(m_data->CvPointsCT.size());
    }

    cv::Point2f center;
    float radius = 0;
    minEnclosingCircle(transformPoints, center, radius);

    cv::Mat centerTransform(3, 1, CV_64F);
    centerTransform.at<double>(0, 0) = center.x;
    centerTransform.at<double>(1, 0) = center.y;
    centerTransform.at<double>(2, 0) = z;

    cv::Mat myCenterMat = transformZ.inv() * centerTransform;
    cv::Point3d myCenter = cv::Point3d(myCenterMat);

    cv::Point3d upPoint = myCenter + radius * kneeHipVector;
    cv::Point3d downPoint = myCenter - radius * kneeHipVector;
    cv::Point3d latPoint = myCenter + radius * medLatVector;
    cv::Point3d medPoint = myCenter - radius * medLatVector;

    cv::Point3d tempPoint;

    circlePoints.push_back(upPoint);

    tempPoint = GetPointOnCircle(upPoint, medPoint, myCenter, 0.2, radius);
    circlePoints.push_back(tempPoint);

    tempPoint = GetPointOnCircle(upPoint, medPoint, myCenter, 0.4, radius);
    circlePoints.push_back(tempPoint);

    tempPoint = GetPointOnCircle(upPoint, medPoint, myCenter, 0.6, radius);
    circlePoints.push_back(tempPoint);

    tempPoint = GetPointOnCircle(upPoint, medPoint, myCenter, 0.8, radius);
    circlePoints.push_back(tempPoint);

    /////////////////////////

    circlePoints.push_back(medPoint);

    tempPoint = GetPointOnCircle(medPoint, downPoint, myCenter, 0.2, radius);
    circlePoints.push_back(tempPoint);

    tempPoint = GetPointOnCircle(medPoint, downPoint, myCenter, 0.4, radius);
    circlePoints.push_back(tempPoint);

    tempPoint = GetPointOnCircle(medPoint, downPoint, myCenter, 0.6, radius);
    circlePoints.push_back(tempPoint);

    tempPoint = GetPointOnCircle(medPoint, downPoint, myCenter, 0.8, radius);
    circlePoints.push_back(tempPoint);

    circlePoints.push_back(downPoint);

    //////////////////

   /* circlePoints.push_back(downPoint);

    tempPoint = GetPointOnCircle(downPoint, medPoint, myCenter, 0.33, radius);
    circlePoints.push_back(tempPoint);

    tempPoint = GetPointOnCircle(downPoint, medPoint, myCenter, 0.67, radius);
    circlePoints.push_back(tempPoint);


    circlePoints.push_back(medPoint);

    tempPoint = GetPointOnCircle(medPoint, upPoint, myCenter, 0.33, radius);
    circlePoints.push_back(tempPoint);

    tempPoint = GetPointOnCircle(medPoint, upPoint, myCenter, 0.67, radius);
    circlePoints.push_back(tempPoint);*/

    /////////////////////////////

    vtkNew<vtkPolyData> polydata;
    polydata->SetPoints(projPoints);

    vtkNew<vtkKdTreePointLocator> kDTree;
    kDTree->SetDataSet(polydata);
    kDTree->BuildLocator();

    for (int i = 0; i < circlePoints.size(); i++)
    {
        double pnt[3];
        pnt[0] = circlePoints[i].x;
        pnt[1] = circlePoints[i].y;
        pnt[2] = circlePoints[i].z;

        vtkIdType tempId = kDTree->FindClosestPoint(pnt);
        double tempDistance = m_data->GetDistance(projPoints->GetPoint(tempId), pnt);

        cv::Point3d vectorRadius = myCenter - circlePoints[i];
        vectorRadius = vectorRadius / sqrt(vectorRadius.dot(vectorRadius));

        cv::Point3d newPoint = circlePoints[i] + tempDistance * vectorRadius;

        circlePoints[i] = newPoint;
        myLines.push_back(RLine(fitPlane.getNormalVector(), newPoint));
    }

    it1 = m_data->CvPointsCT.begin();
    it2 = m_data->CvPointsCT.end();

    std::vector<double> distances(circlePoints.size(), 9999999.0);

    for (; it1 != it2; ++it1)
    {
        for (int i = 0; i < circlePoints.size(); i++)
        {
            double dist = myLines[i].getDistanceFromPoint(*it1);
            if (dist < distances[i])
            {
                distances[i] = dist;
                circlePoints[i] = *it1;
            }
        }
    }

    for (int i = 0; i < circlePoints.size(); i++)
    {
        mPoints.push_back(m_data->cvPointToItk(circlePoints[i]));
    }

}

std::vector<PointTypeITK> KneeCapRegistration::GetRegistrationPoints() const
{
    return mPoints;
}

/*
bool KneeCapRegistration::RegistrationLandmarks(const PointTypeITK& pTopPointCamera, const PointTypeITK& pDownPointCamera, const PointTypeITK& pLateralPointCamera, const PointTypeITK& pMedialPointCamera, double& error)
{
    if (isVTK == true)
    {
        return RegistrationLandmarksVTK(pTopPointCamera, pDownPointCamera, pLateralPointCamera, pMedialPointCamera, error);
    }
    else
    {
        return RegistrationLandmarksPCL(pTopPointCamera, pDownPointCamera, pLateralPointCamera, pMedialPointCamera, error);
    }
}
*/

bool KneeCapRegistration::MakeRegistration(const std::vector<PointTypeITK>& pBonePoints)
{
    return MakeRegistrationLS(pBonePoints);

    /*if (useRandomAlignment == true)
    {
        return MakeRegistrationLSRandom(pBonePoints, pTopPointCamera, pDownPointCamera, pLateralPointCamera, pMedialPointCamera);
    }
    else
    {
        return MakeRegistrationLS(pBonePoints, pTopPointCamera, pDownPointCamera, pLateralPointCamera, pMedialPointCamera);
    }*/
}

/*
bool KneeCapRegistration::RegistrationLandmarksPCL(const PointTypeITK& pTopPointCamera, const PointTypeITK& pDownPointCamera, const PointTypeITK& pLateralPointCamera, const PointTypeITK& pMedialPointCamera, double& error)
{
    double near_radius = 5.0;

    std::vector<PointTypeITK> source, target;
    PointTypeITK top_temp_ct, down_temp_ct, lateral_temp_ct, medial_temp_ct;

    std::vector<pcl::PointXYZ> top_list_ct;
    std::vector<pcl::PointXYZ> down_list_ct;
    std::vector<pcl::PointXYZ> lateral_list_ct;
    std::vector<pcl::PointXYZ> medial_list_ct;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_all;
    pcl::PointCloud<pcl::PointXYZ>::Ptr landMarkPoints(new pcl::PointCloud<pcl::PointXYZ>);

    std::vector<pcl::PointXYZ> pclPoints;

    pclPoints.push_back(m_data->itkPointToPCL(topPointCT));
    pclPoints.push_back(m_data->itkPointToPCL(downPointCT));
    pclPoints.push_back(m_data->itkPointToPCL(lateralPointCT));
    pclPoints.push_back(m_data->itkPointToPCL(medialPointCT));

    cloud_all = m_data->getNearPointsToReferencePoints(pclPoints, top_list_ct, down_list_ct, lateral_list_ct, medial_list_ct, near_radius);

    if (top_list_ct.size() > 0 && down_list_ct.size() > 0 && lateral_list_ct.size() > 0 && medial_list_ct.size() > 0)
    {
        top_temp_ct = m_data->pclPointToItk(top_list_ct[0]);
        down_temp_ct = m_data->pclPointToItk(down_list_ct[0]);
        lateral_temp_ct = m_data->pclPointToItk(lateral_list_ct[0]);
        medial_temp_ct = m_data->pclPointToItk(medial_list_ct[0]);
    }
    else
    {
        throw RegistrationException("Check your CT landmarks. Distance from the bone model is too large.");
        return false;
    }

    source.push_back(pTopPointCamera);
    source.push_back(pDownPointCamera);
    source.push_back(pLateralPointCamera);
    source.push_back(pMedialPointCamera);

    target.push_back(top_temp_ct);
    target.push_back(down_temp_ct);
    target.push_back(lateral_temp_ct);
    target.push_back(medial_temp_ct);

    Eigen::Matrix4d rigidTransform = m_data->InitRegistration(source, target);

    std::vector<PointTypeITK> cameraPoints = { pTopPointCamera, pDownPointCamera, pLateralPointCamera, pMedialPointCamera };

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


bool KneeCapRegistration::RegistrationLandmarksVTK(const PointTypeITK& pTopPointCamera, const PointTypeITK& pDownPointCamera, const PointTypeITK& pLateralPointCamera, const PointTypeITK& pMedialPointCamera, double& error)
{
    vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    implicitPolyDataDistance->SetInput(poly);

    std::vector<PointTypeITK> source, target;
    PointTypeITK my_top_ct, my_down_ct, my_lateral_ct, my_medial_ct;

    double signedDistance1 = m_data->GetNearPoint(implicitPolyDataDistance, topPointCT, my_top_ct);
    double signedDistance2 = m_data->GetNearPoint(implicitPolyDataDistance, downPointCT, my_down_ct);
    double signedDistance3 = m_data->GetNearPoint(implicitPolyDataDistance, lateralPointCT, my_lateral_ct);
    double signedDistance4 = m_data->GetNearPoint(implicitPolyDataDistance, medialPointCT, my_medial_ct);

    if (signedDistance1 > 5 || signedDistance2 > 5 || signedDistance3 > 5 || signedDistance4 > 5)
    {
        throw RegistrationException("Check your CT landmarks. Distance from the bone model is too large.");
        return false;
    }

    source.push_back(pTopPointCamera);
    source.push_back(pDownPointCamera);
    source.push_back(pLateralPointCamera);
    source.push_back(pMedialPointCamera);

    target.push_back(my_top_ct);
    target.push_back(my_down_ct);
    target.push_back(my_lateral_ct);
    target.push_back(my_medial_ct);

    Eigen::Matrix4d rigidTransform = m_data->InitRegistration(source, target);

    PointTypeITK temp;

    double a = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pTopPointCamera, temp);
    double b = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pDownPointCamera, temp);
    double c = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pLateralPointCamera, temp);
    double d = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pMedialPointCamera, temp);

    std::vector<double> distances = { a, b, c, d};
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
*/
bool KneeCapRegistration::MakeRegistrationLS(const std::vector<PointTypeITK>& pBonePoints)
{
    cv::Mat data = Registration::GetTranslationRotation(pBonePoints, mPoints);

    LeastSquaresICP myICP(pBonePoints);

    double error;

    /*if (isVTK == true)
    {
        error = myICP.LeastSquares(Registration::poly, data);
    }
    else
    {
        error = myICP.LeastSquares(m_data->PointsCT, data);
    }
*/
	vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
	implicitPolyDataDistance->SetInput(Registration::poly);
    error = myICP.LeastSquares(implicitPolyDataDistance, data);

    Registration::MakeResult(data, error);

    return true;
}

/*
bool KneeCapRegistration::MakeRegistrationLSRandom(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pTopPointCamera, const PointTypeITK& pDownPointCamera, const PointTypeITK& pLateralPointCamera, const PointTypeITK& pMedialPointCamera)
{
    double tolerance = 10.0;
    double near_radius = 5.0;
    double max_distance = tolerance;

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_all(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr bone_points(new pcl::PointCloud<pcl::PointXYZ>);

    (*cloud_all) = *(m_data->PointsCT);

    if (pBonePoints.size() == 0)
    {
        return false;
    }
    else
    {
        m_data->itkPointVectorToPCLPointVector(pBonePoints, bone_points);
    }

    std::vector<pcl::PointXYZ> top_list_ct;
    std::vector<pcl::PointXYZ> down_list_ct;
    std::vector<pcl::PointXYZ> lateral_list_ct;
    std::vector<pcl::PointXYZ> medial_list_ct;

    std::vector<pcl::PointXYZ> pclPoints;

    pclPoints.push_back(m_data->itkPointToPCL(topPointCT));
    pclPoints.push_back(m_data->itkPointToPCL(downPointCT));
    pclPoints.push_back(m_data->itkPointToPCL(lateralPointCT));
    pclPoints.push_back(m_data->itkPointToPCL(medialPointCT));

    m_data->getMostNearPointsToReferencePoints(pclPoints, top_list_ct, down_list_ct, lateral_list_ct, medial_list_ct, near_radius);

    Eigen::Matrix4d rigidTransform;

    std::vector<PointTypeITK> source, target;

    double radius = 2.0;

    PointTypeITK top_temp_ct, down_temp_ct, lateral_temp_ct, medial_temp_ct;

    if (top_list_ct.size() > 0 && down_list_ct.size() > 0 && lateral_list_ct.size() > 0 && medial_list_ct.size() > 0)
    {
        top_temp_ct = m_data->pclPointToItk(top_list_ct[0]);
        down_temp_ct = m_data->pclPointToItk(down_list_ct[0]);
        lateral_temp_ct = m_data->pclPointToItk(lateral_list_ct[0]);
        medial_temp_ct = m_data->pclPointToItk(medial_list_ct[0]);
    }
    else
    {
        throw RegistrationException("Check your CT landmarks. Distance from the bone model is too large.");
        return false;
    }

    PointTypeITK my_top_ct = top_temp_ct;
    PointTypeITK my_down_ct = down_temp_ct;
    PointTypeITK my_lateral_ct = lateral_temp_ct;
    PointTypeITK my_medial_ct = medial_temp_ct;

    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(cloud_all);

    for (int i = 0; i < 10000; i++)
    {
        source.clear();
        target.clear();

        source.push_back(pTopPointCamera);
        source.push_back(pDownPointCamera);
        source.push_back(pLateralPointCamera);
        source.push_back(pMedialPointCamera);

        target.push_back(my_top_ct);
        target.push_back(my_down_ct);
        target.push_back(my_lateral_ct);
        target.push_back(my_medial_ct);

        rigidTransform = m_data->InitRegistration(source, target);

        pcl::PointCloud<pcl::PointXYZ>::Ptr output_cloud(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::transformPointCloud(*bone_points, *output_cloud, rigidTransform);

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
            top_temp_ct = my_top_ct;
            down_temp_ct = my_down_ct;
            lateral_temp_ct = my_lateral_ct;
            medial_temp_ct = my_medial_ct;
        }

        my_top_ct = m_data->getNearPointRandom(kdtree, cloud_all, getPointInsideSphere(top_temp_ct, radius));
        my_down_ct = m_data->getNearPointRandom(kdtree, cloud_all, getPointInsideSphere(down_temp_ct, radius));
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

    LeastSquaresICP myICP(pBonePoints);

    if (myICP.getMaxError() > m_data->mError)
    {
        myICP.setMaxError(m_data->mError);
    }

    double error;

    //std::cout << "My Random Error: " << m_data->mError << std::endl;

    if (isVTK == true)
    {
        error = myICP.LeastSquares(poly, data);
    }
    else
    {
        error = myICP.LeastSquares(m_data->PointsCT, data);
    }

    //std::cout << "My ICP Error: " << error << std::endl;

    if (m_data->mError > error)
    {
        Registration::MakeResult(data, error);
    }

    return true;
}
*/

cv::Point3d KneeCapRegistration::GetPointOnCircle(const cv::Point3d& a, const cv::Point3d& b, const cv::Point3d& center, double perCent, double radius) const
{
    cv::Point3d vector = (a + perCent * (b - a)) - center;
    vector = vector / sqrt(vector.dot(vector));
    return (center + radius * vector);
}

cv::Mat KneeCapRegistration::GetRotateZ(const cv::Point3d& vector) const
{
    cv::Mat vectorMat(3, 1, CV_64F);
    vectorMat.at<double>(0, 0) = vector.x;
    vectorMat.at<double>(1, 0) = vector.y;
    vectorMat.at<double>(2, 0) = vector.z;

    cv::Point3d normalXY(0, 0, 1);
    cv::Point3d rotationAxis = normalXY.cross(vector);
    rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
    double rotationAngle = RLine::getAngleBetweenVectors(normalXY, vector);
    cv::Mat rotation_1 = RLine::getRotateMatrix(rotationAxis, -rotationAngle);
    cv::Mat rotation_2 = RLine::getRotateMatrix(rotationAxis, rotationAngle);
    cv::Mat rotateVector_1 = rotation_1 * vectorMat;
    cv::Mat rotateVector_2 = rotation_2 * vectorMat;

    cv::Mat rotate;

    double distance_1 = RLine::getAngleBetweenVectors(cv::Point3d(rotateVector_1), normalXY);
    double distance_2 = RLine::getAngleBetweenVectors(cv::Point3d(rotateVector_2), normalXY);

    if (distance_1 < distance_2)
    {
        rotate = rotation_1;
    }
    else
    {
        rotate = rotation_2;
    }

    return rotate;
}
