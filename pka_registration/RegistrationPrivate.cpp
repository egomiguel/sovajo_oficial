#include "RegistrationPrivate.hpp"
#include <random>
//#include <pcl/recognition/ransac_based/trimmed_icp.h>
#include "RPlane.hpp"
#include "LeastSquaresICP.hpp"

using namespace PKA::REGISTRATION;

RegistrationPrivate::RegistrationPrivate()
{
    //PointsCT = pcl::PointCloud<pcl::PointXYZ>::Ptr(new pcl::PointCloud<pcl::PointXYZ>);
    mMatrix = Eigen::Matrix4d::Identity();
    mError = 9999999;
}


Eigen::Matrix4d RegistrationPrivate::InitRegistration(const std::vector<PointTypeITK>& source, const std::vector<PointTypeITK>& target)
{
    Eigen::Matrix4d transform;
    
    if (source.size() != target.size() || source.size() == 0)
    {
        return Eigen::Matrix4d::Identity();
    }

    int tSize = source.size();

    ////////////////////////////////////////////////////////////// Find Centroid
    Eigen::Matrix<double, 4, 1> centroid_src, centroid_tgt;

    centroid_src.setZero();
    centroid_tgt.setZero();

    auto it1 = source.begin();
    auto it2 = source.end();

    for (; it1 != it2; ++it1)
    {
        centroid_src[0] += (*it1)[0];
        centroid_src[1] += (*it1)[1];
        centroid_src[2] += (*it1)[2];
    }

    centroid_src /= static_cast<double>(tSize);
    centroid_src[3] = 1;

    it1 = target.begin();
    it2 = target.end();

    for (; it1 != it2; ++it1)
    {
        centroid_tgt[0] += (*it1)[0];
        centroid_tgt[1] += (*it1)[1];
        centroid_tgt[2] += (*it1)[2];
    }

    centroid_tgt /= static_cast<double>(tSize);
    centroid_tgt[3] = 1;
    
    //////////////////////////////////////////////////////////////////////////

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> cloud_src_demean, cloud_tgt_demean;
    cloud_src_demean = Eigen::Matrix<double, 4, Eigen::Dynamic>::Zero(4, tSize);
    cloud_tgt_demean = Eigen::Matrix<double, 4, Eigen::Dynamic>::Zero(4, tSize);

    it1 = source.begin();
    it2 = source.end();
    int i = 0;
    while (it1 != it2)
    {
        cloud_src_demean(0, i) = (*it1)[0] - centroid_src[0];
        cloud_src_demean(1, i) = (*it1)[1] - centroid_src[1];
        cloud_src_demean(2, i) = (*it1)[2] - centroid_src[2];
        ++i;
        ++it1;
    }

    it1 = target.begin();
    it2 = target.end();
    i = 0;
    while (it1 != it2)
    {
        cloud_tgt_demean(0, i) = (*it1)[0] - centroid_tgt[0];
        cloud_tgt_demean(1, i) = (*it1)[1] - centroid_tgt[1];
        cloud_tgt_demean(2, i) = (*it1)[2] - centroid_tgt[2];
        ++i;
        ++it1;
    }

    ////////////////////////////////////////////////////////////////////

    transform = Eigen::Matrix4d::Identity();

    // Assemble the correlation matrix H = source * target'
    Eigen::Matrix<double, 3, 3> H = (cloud_src_demean * cloud_tgt_demean.transpose()).topLeftCorner(3, 3);

    // Compute the Singular Value Decomposition
    Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3> > svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix<double, 3, 3> u = svd.matrixU();
    Eigen::Matrix<double, 3, 3> v = svd.matrixV();

    // Compute R = V * U'
    if (u.determinant() * v.determinant() < 0)
    {
        for (int x = 0; x < 3; ++x)
            v(x, 2) *= -1;
    }

    Eigen::Matrix<double, 3, 3> R = v * u.transpose();

    // Return the correct transformation
    transform.topLeftCorner(3, 3) = R;
    const Eigen::Matrix<double, 3, 1> Rc(R * centroid_src.head(3));
    transform.block(0, 3, 3, 1) = centroid_tgt.head(3) - Rc;

    //////////////////////////////////////////////////////////////////////

    /*
    pcl::PointCloud<pcl::PointXYZ> cloud_in;
    pcl::PointCloud<pcl::PointXYZ> cloud_out;
    pcl::registration::TransformationEstimationSVD<pcl::PointXYZ, pcl::PointXYZ, double> TESVD;

    for (int i = 0; i < tSize; i++)
    {
        cloud_in.points.push_back(itkPointToPCL(source[i]));
        cloud_out.points.push_back(itkPointToPCL(target[i]));
    }

    TESVD.estimateRigidTransformation(cloud_in, cloud_out, transform);*/

    return transform;
}


double RegistrationPrivate::GetNearPoint(const vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance, const PointTypeITK& point, PointTypeITK& nearPoint)
{
    double myClosest[3];
    double pnt[3];
    pnt[0] = point[0];
    pnt[1] = point[1];
    pnt[2] = point[2];

    double distance = implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);

    nearPoint[0] = myClosest[0];
    nearPoint[1] = myClosest[1];
    nearPoint[2] = myClosest[2];

    return abs(distance);
}

void RegistrationPrivate::GetKneeCapPoints(const vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance, const cv::Point3d& a, const cv::Point3d& b, const cv::Point3d& vector, std::vector<PointTypeITK>& pPoints)
{
    double k = 50;

    cv::Point3d p1 = (a + 0.25 * (b - a)) + k * vector;
    cv::Point3d p2 = (a + 0.75 * (b - a)) + k * vector;

    double myClosest[3];
    double pnt[3];

    pnt[0] = p1.x;
    pnt[1] = p1.y;
    pnt[2] = p1.z;

    implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);

    pPoints.push_back(vtkPointToItk(myClosest));

    pnt[0] = p2.x;
    pnt[1] = p2.y;
    pnt[2] = p2.z;

    implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);

    pPoints.push_back(vtkPointToItk(myClosest));
}

double RegistrationPrivate::TransformPointDistance(const vtkSmartPointer<vtkImplicitPolyDataDistance> implicitDistance, const Eigen::Matrix4d& matrix, const PointTypeITK& point, PointTypeITK& pointOut)
{
    Eigen::Vector4d myPoint;
    myPoint << point[0], point[1], point[2], 1;
    Eigen::Vector4d myTransform = matrix * myPoint;
    double result[3];
    result[0] = myTransform(0);
    result[1] = myTransform(1);
    result[2] = myTransform(2);

    double myClosest[3];

    double distance = implicitDistance->EvaluateFunctionAndGetClosestPoint(result, myClosest);

    pointOut[0] = myTransform[0];
    pointOut[1] = myTransform[1];
    pointOut[2] = myTransform[2];

    return abs(distance);
}

double RegistrationPrivate::GetITKPointDistance(const PointTypeITK& point1, const PointTypeITK& point2)
{
    cv::Point3d diff = { point1[0] - point2[0], point1[1] - point2[1] , point1[2] - point2[2] };
    return sqrt(diff.dot(diff));
}

/*
void RegistrationPrivate::alignCloud(const pcl::PointCloud<pcl::PointXYZ>::Ptr source_cloud, const pcl::PointCloud<pcl::PointXYZ>::Ptr target_cloud, double near_radius)
{
    double radius = floor(near_radius) + 1.0;
    double max_angle = radius * 2.0;
    double max_translation = radius * 10.0;
    int set_size = (max_translation * 4) + 1;
    double compare_dis = near_radius;


    std::vector<double> translation;
    translation.reserve(set_size);
    double trasnlation_step = (2.0 * max_translation) / double(set_size - 1);
    for (double i = -max_translation; i < max_translation + trasnlation_step; i += trasnlation_step)
    {
        translation.push_back(i);
        if (translation.size() == set_size)
        {
            break;
        }
    }

    std::vector<double> rotation;
    rotation.reserve(set_size);
    double rotation_step = (2.0 * max_angle) / double(set_size - 1);
    for (double i = -max_angle; i < max_angle + rotation_step; i += rotation_step)
    {
        rotation.push_back(i*M_PI / 180.0);
        if (rotation.size() == set_size)
        {
            break;
        }
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> matrix_distr(0, set_size - 1);
    int thx_pos, thy_pos, thz_pos, x_pos, y_pos, z_pos;
    double thX, thY, thZ, x, y, z;

    auto ICP = pcl::recognition::TrimmedICP<pcl::PointXYZ, double>::TrimmedICP();
    //Eigen::Matrix4f guess = Eigen::Matrix4f::Identity();
    ICP.init(target_cloud);

    Eigen::Matrix4d align = Eigen::Matrix4d::Identity();
    int source_cloud_size = (*source_cloud).size();

    for (int i = 0; i < 1000; i++)
    {
        thx_pos = matrix_distr(gen);
        thy_pos = matrix_distr(gen);
        thz_pos = matrix_distr(gen);

        thX = rotation[thx_pos];
        thY = rotation[thy_pos];
        thZ = rotation[thz_pos];

        x_pos = matrix_distr(gen);
        y_pos = matrix_distr(gen);
        z_pos = matrix_distr(gen);

        x = translation[x_pos];
        y = translation[y_pos];
        z = translation[z_pos];

        //std::cout << x << " " << y << " " << z << " " << thX << " " << thY << " " << thZ << std::endl;

        Eigen::AngleAxisd init_rotationZ(thZ, Eigen::Vector3d::UnitZ());
        Eigen::AngleAxisd init_rotationY(thY, Eigen::Vector3d::UnitY());
        Eigen::AngleAxisd init_rotationX(thX, Eigen::Vector3d::UnitX());

        Eigen::Translation3d init_translation(x, y, z);

        Eigen::Matrix4d init_guess = (init_translation * init_rotationX * init_rotationY * init_rotationZ).matrix();

        ICP.align(*source_cloud, source_cloud_size - 1, init_guess);

        pcl::PointCloud<pcl::PointXYZ>::Ptr output_cloud(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::transformPointCloud(*source_cloud, *output_cloud, init_guess);

        double dis = getFitnessScore_max(target_cloud, output_cloud);

        if (dis < compare_dis)
        {
            align = init_guess;
            compare_dis = dis;
            std::cout << "Align error ICP: " << compare_dis << std::endl;
            if (compare_dis <= 0.1)
            {
                //std::cout << "Second number Iteration: " << i << std::endl;
                break;
            }
        }
    }

    mMatrix = align * mMatrix;
    mError = compare_dis;
    ///////////////////////////////////////////////////////
    
    //Eigen::Matrix3d ab;
    //ab(0, 0) = mMatrix(0, 0);
    //ab(0, 1) = mMatrix(0, 1);
    //ab(0, 2) = mMatrix(0, 2);

    //ab(1, 0) = mMatrix(1, 0);
    //ab(1, 1) = mMatrix(1, 1);
    //ab(1, 2) = mMatrix(1, 2);

    //ab(2, 0) = mMatrix(2, 0);
    //ab(2, 1) = mMatrix(2, 1);
    //ab(2, 2) = mMatrix(2, 2);

    //Eigen::Vector3d teuler = ab.eulerAngles(0, 1, 2);
    //teuler = (teuler / M_PI)*180.0;
    //cout << "Angles rotation around x y z axis:\n" << teuler << endl;
    //Eigen::Vector3d tran;
    //tran << mMatrix(0, 3), mMatrix(1, 3), mMatrix(2, 3);
    //cout << "Translation x y z:\n" << tran << endl;
    ////////////////////////////////////////////////////////////
}
*/

/*
double RegistrationPrivate::getFitnessScore_max(const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_sc,
    const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_targe, bool show)
{
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(cloud_sc);

    int K = 1;

    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);
    double maxDis = 0.0;
    for (size_t i = 0; i < cloud_targe->points.size(); ++i)
    {
        if (kdtree.nearestKSearch(cloud_targe->points[i], K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
        {
            if (show == true)
            {
                std::cout << cloud_sc->points[pointIdxNKNSearch[0]] << " " << cloud_targe->points[i] << " " << sqrt(pointNKNSquaredDistance[0]) << std::endl;
            }

            if (pointNKNSquaredDistance[0] > maxDis)
            {
                maxDis = pointNKNSquaredDistance[0];
            }
        }
    }

    return sqrt(maxDis);
}
*/

/*
double RegistrationPrivate::getFitnessScore_max(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_targe)
{
    int K = 1;

    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);
    double maxDis = 0.0;
    for (size_t i = 0; i < cloud_targe->points.size(); ++i)
    {
        if (kdtree.nearestKSearch(cloud_targe->points[i], K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
        {
            if (pointNKNSquaredDistance[0] > maxDis)
            {
                maxDis = pointNKNSquaredDistance[0];
            }
        }
        pointIdxNKNSearch.clear();
        pointNKNSquaredDistance.clear();
    }

    return sqrt(maxDis);
}
*/

double RegistrationPrivate::getFitnessScore_max(const vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance, const std::vector<Eigen::Vector4d>& cloud_targe)
{
    auto it1 = cloud_targe.begin();
    auto it2 = cloud_targe.end();
    double distance = -1.0;

    for (; it1 != it2; ++it1)
    {
        double myClosest[3];
        double pnt[3];
        pnt[0] = (*it1)[0];
        pnt[1] = (*it1)[1];
        pnt[2] = (*it1)[2];

        double distanceTemp = implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);

        if (abs(distanceTemp) > distance)
        {
            distance = abs(distanceTemp);
        }
    }

    return distance;
}

/*
void RegistrationPrivate::reduceBiggerCloud(const pcl::PointCloud<pcl::PointXYZ>::Ptr biggerCloud, const pcl::PointCloud<pcl::PointXYZ>::Ptr smallerCloud,
    pcl::PointCloud<pcl::PointXYZ>::Ptr reduceBiggerCloudOut, double radius)
{
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(biggerCloud);


    std::vector<int> pointsIdx;
    std::vector<float> pointSquaredDistance;
    std::set<int> points;

    auto it1 = smallerCloud->points.begin();
    auto it2 = smallerCloud->points.end();

    std::vector<int>::iterator itPointBegin, itPointEnd;

    for (; it1 != it2; ++it1)
    {

        if (kdtree.radiusSearch((*it1), radius, pointsIdx, pointSquaredDistance))
        {
            itPointBegin = pointsIdx.begin();
            itPointEnd = pointsIdx.end();

            for (; itPointBegin != itPointEnd; ++itPointBegin)
            {
                points.insert(*itPointBegin);
            }
            pointsIdx.clear();
            pointSquaredDistance.clear();
        }
    }

    auto setIt1 = points.begin();
    auto setIt2 = points.end();

    for (; setIt1 != setIt2; ++setIt1)
    {
        reduceBiggerCloudOut->points.push_back(biggerCloud->points[*setIt1]);
    }
}
*/

/*
void RegistrationPrivate::reduceBiggerCloud(const pcl::PointCloud<pcl::PointXYZ>::Ptr biggerCloud, const pcl::PointCloud<pcl::PointXYZ>::Ptr smallerCloud,
    Eigen::MatrixXd& outPut, double radius)
{
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(biggerCloud);


    std::vector<int> pointsIdx;
    std::vector<float> pointSquaredDistance;
    std::set<int> points;

    auto it1 = smallerCloud->points.begin();
    auto it2 = smallerCloud->points.end();

    std::vector<int>::iterator itPointBegin, itPointEnd;

    for (; it1 != it2; ++it1)
    {

        if (kdtree.radiusSearch((*it1), radius, pointsIdx, pointSquaredDistance))
        {
            itPointBegin = pointsIdx.begin();
            itPointEnd = pointsIdx.end();

            for (; itPointBegin != itPointEnd; ++itPointBegin)
            {
                points.insert(*itPointBegin);
            }
            pointsIdx.clear();
            pointSquaredDistance.clear();
        }
    }

    auto setIt1 = points.begin();
    auto setIt2 = points.end();

    Eigen::MatrixXd outPutTemp(points.size(), 3);

    int cont = 0;
    for (; setIt1 != setIt2; ++setIt1)
    {
        pcl::PointXYZ myPoint = biggerCloud->points[*setIt1];
        outPutTemp(cont, 0) = myPoint.x;
        outPutTemp(cont, 1) = myPoint.y;
        outPutTemp(cont, 2) = myPoint.z;
        cont++;
    }
    outPut = outPutTemp;
}
*/

Eigen::Vector3d RegistrationPrivate::itkPointToEigen(const PointTypeITK& point)
{
    Eigen::Vector3d vec;
    vec << point[0], point[1], point[2];
    return vec;
}


/*
void RegistrationPrivate::itkPointVectorToPCLPointVector(const std::vector<PointTypeITK>& pointsIn, pcl::PointCloud<pcl::PointXYZ>::Ptr pointsOut)
{
    auto it1 = pointsIn.begin();
    auto it2 = pointsIn.end();
    pcl::PointXYZ pclPoint;
    for (; it1 != it2; ++it1)
    {
        pointsOut->points.push_back(itkPointToPCL(*it1));
    }
}
*/

void RegistrationPrivate::itkPointVectorToEigenPointVector4(const std::vector<PointTypeITK>& pointsIn, std::vector<Eigen::Vector4d>& pointsOut)
{
    auto it1 = pointsIn.begin();
    auto it2 = pointsIn.end();

    for (; it1 != it2; ++it1) 
    {
        Eigen::Vector4d point4;
        point4 << (*it1)[0], (*it1)[1], (*it1)[2], 1;
        pointsOut.push_back(point4);
    }
}

void RegistrationPrivate::transformVectorEigen4(const std::vector<Eigen::Vector4d>& pointsIn, std::vector<Eigen::Vector4d>& pointsOut, const Eigen::Matrix4d& pMatrix)
{
    auto it1 = pointsIn.begin();
    auto it2 = pointsIn.end();

    for (; it1 != it2; ++it1)
    {
        Eigen::Vector4d point4 = pMatrix * (*it1);
        pointsOut.push_back(point4);
    }
}

cv::Point3d RegistrationPrivate::itkPointToCV(const PointTypeITK& point)
{
    cv::Point3d cvPoint;
    cvPoint.x = point[0];
    cvPoint.y = point[1];
    cvPoint.z = point[2];
    return cvPoint;
}
/*
PointTypeITK RegistrationPrivate::pclPointToItk(const pcl::PointXYZ& point)
{
    PointTypeITK itkPoint;
    itkPoint[0] = point.x;
    itkPoint[1] = point.y;
    itkPoint[2] = point.z;
    return itkPoint;
}
*/
Eigen::Vector3d RegistrationPrivate::cvPointToEigen(const cv::Point3d& point)
{
    Eigen::Vector3d newPoint;
    newPoint << point.x, point.y, point.z;
    return newPoint;
}

/*
pcl::PointCloud<pcl::PointXYZ>::Ptr RegistrationPrivate::getNearPointsToReferencePoints(const std::vector<pcl::PointXYZ>& points, std::vector<pcl::PointXYZ>& point_list_1, std::vector<pcl::PointXYZ>& point_list_2,
    std::vector<pcl::PointXYZ>& point_list_3, std::vector<pcl::PointXYZ>& point_list_4, double radius)
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_all(new pcl::PointCloud<pcl::PointXYZ>);

    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(PointsCT);

    std::vector<int> pointsIdx;
    std::vector<float> pointSquaredDistance;
    std::vector<int>::iterator it1, it2;
    std::set<int> ids;
    //float LeafSize = 0.1;

    
    for (int i = 0; i < points.size(); i++)
    {
        if (kdtree.radiusSearch(points[i], radius, pointsIdx, pointSquaredDistance))
        {
            it1 = pointsIdx.begin();
            it2 = pointsIdx.end();

            for (; it1 != it2; ++it1)
            {
                ids.insert(*it1);
            }
            pointsIdx.clear();
            pointSquaredDistance.clear();

            if (kdtree.nearestKSearch(points[i], 1, pointsIdx, pointSquaredDistance))
            {
                it1 = pointsIdx.begin();
                
                if (i == 0)
                {
                    point_list_1.push_back(PointsCT->points[*it1]);
                }
                else if (i == 1)
                {
                    point_list_2.push_back(PointsCT->points[*it1]);
                }
                else if (i == 2)
                {
                    point_list_3.push_back(PointsCT->points[*it1]);
                }
                else if (i == 3)
                {
                    point_list_4.push_back(PointsCT->points[*it1]);
                }
            }

            pointsIdx.clear();
            pointSquaredDistance.clear();

            //pcl::VoxelGrid<pcl::PointXYZ> filter;
            //filter.setInputCloud(point_list_1_temp);
            //filter.setLeafSize(LeafSize, LeafSize, LeafSize);
            //filter.filter(*point_list_1);
        }
    }
    
    std::set<int>::iterator itSet1 = ids.begin();
    std::set<int>::iterator itSet2 = ids.end();

    for ( ; itSet1 != itSet2; ++itSet1 )
    {
        cloud_all->points.push_back(PointsCT->points[*itSet1]);
    }

    return cloud_all;
}
*/

/*
PointTypeITK RegistrationPrivate::getNearPointRandom(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, const pcl::PointCloud<pcl::PointXYZ>::Ptr pointList, PointTypeITK& point)
{
    std::vector<int> pointsIdx;
    std::vector<float> pointSquaredDistance;
    std::vector<int>::iterator it1;
    int K = 1;

    pcl::PointXYZ pointPCL, nearPointPCL;
    pointPCL.x = point[0];
    pointPCL.y = point[1];
    pointPCL.z = point[2];

    if (kdtree.nearestKSearch(pointPCL, K, pointsIdx, pointSquaredDistance))
    {
        it1 = pointsIdx.begin();

        nearPointPCL = pointList->points[*it1];
    }
    PointTypeITK result;
    result[0] = nearPointPCL.x;
    result[1] = nearPointPCL.y;
    result[2] = nearPointPCL.z;
    return result;
}
*/

/*
PointTypeITK RegistrationPrivate::getNearPointRandom(const pcl::KdTreeFLANN<pcl::PointXYZ>& kdtree, PointTypeITK& point, double radius)
{
    std::vector<int> pointsIdx;
    std::vector<float> pointSquaredDistance;

    pcl::PointXYZ pointPCL;
    pointPCL.x = point[0];
    pointPCL.y = point[1];
    pointPCL.z = point[2];

    if (kdtree.radiusSearch(pointPCL, radius, pointsIdx, pointSquaredDistance))
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distr(pointsIdx.size() / 3, pointsIdx.size() - 1);
        pcl::PointXYZ nearPoint = PointsCT->points[distr(gen)];

        PointTypeITK myPoint;
        myPoint[0] = nearPoint.x;
        myPoint[1] = nearPoint.y;
        myPoint[2] = nearPoint.z;

        return myPoint;
    }
    else
    {
        return point;
    }

}
*/

/*
void RegistrationPrivate::getMostNearPointsToReferencePoints(const std::vector<pcl::PointXYZ>& points, std::vector<pcl::PointXYZ>& point_list_1, std::vector<pcl::PointXYZ>& point_list_2,
    std::vector<pcl::PointXYZ>& point_list_3, std::vector<pcl::PointXYZ>& point_list_4, double radius)
{
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(PointsCT);

    std::vector<int> pointsIdx;
    std::vector<float> pointSquaredDistance;
    std::vector<int>::iterator it1;
    int K = 1;
    double squareRadius = radius * radius;

    for (int i = 0; i < points.size(); i++)
    {
        if (kdtree.nearestKSearch(points[i], K, pointsIdx, pointSquaredDistance))
        {
            it1 = pointsIdx.begin();
            if (pointSquaredDistance[0] < squareRadius)
            {
                if (i == 0)
                {
                    point_list_1.push_back(PointsCT->points[*it1]);
                }
                else if (i == 1)
                {
                    point_list_2.push_back(PointsCT->points[*it1]);
                }
                else if (i == 2)
                {
                    point_list_3.push_back(PointsCT->points[*it1]);
                }
                else if (i == 3)
                {
                    point_list_4.push_back(PointsCT->points[*it1]);
                }
            }
            pointsIdx.clear();
            pointSquaredDistance.clear();
        }
    }
}
*/
/*
PointTypeITK RegistrationPrivate::getRandomPoint(const pcl::PointCloud<pcl::PointXYZ>::Ptr pointList, const PointTypeITK& point)
{
    PointTypeITK result;
    pcl::PointXYZ result_temp;
    pcl::PointXYZ newPoint;
    newPoint.x = point[0];
    newPoint.y = point[1];
    newPoint.z = point[2];

    std::random_device rd;
    std::mt19937 gen(rd());

    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(pointList);

    std::vector<int> pointsIdx;
    std::vector<float> pointSquaredDistance;
    int pos;

    if (kdtree.radiusSearch(newPoint, 2.0, pointsIdx, pointSquaredDistance))
    {
        std::uniform_int_distribution<> distr(0, pointsIdx.size() - 1);
        pos = distr(gen);
        result_temp = pointList->points[pointsIdx[pos]];
    }
    result[0] = result_temp.x;
    result[1] = result_temp.y;
    result[2] = result_temp.z;

    return result;
}
*/

cv::Mat RegistrationPrivate::cvPointToMat(const cv::Point3d& point)
{
    cv::Mat mat(3, 1, CV_64F);
    mat.at<double>(0, 0) = point.x;
    mat.at<double>(1, 0) = point.y;
    mat.at<double>(2, 0) = point.z;
    return mat;
}


PointTypeITK RegistrationPrivate::cvPointToItk(const cv::Point3d& point)
{
    PointTypeITK pointITK;
    pointITK[0] = point.x;
    pointITK[1] = point.y;
    pointITK[2] = point.z;
    return pointITK;
}

PointTypeITK RegistrationPrivate::vtkPointToItk(const double point[3])
{
    PointTypeITK pointITK;
    pointITK[0] = point[0];
    pointITK[1] = point[1];
    pointITK[2] = point[2];
    return pointITK;
}

double RegistrationPrivate::GetDistance(const double p1[3], const double p2[3]) const
{
    double a = abs(p1[0] - p2[0]);
    double b = abs(p1[1] - p2[1]);
    double c = abs(p1[2] - p2[2]);

    return sqrt(a * a + b * b + c * c);
}

double RegistrationPrivate::GetDistance(const cv::Point3d& p1, const cv::Point3d& p2) const
{
    cv::Point3d diff = p2 - p1;
    return sqrt(diff.dot(diff));
}

cv::Mat RegistrationPrivate::GetRotateZ(const cv::Point3d& vector)
{
    cv::Point3d normalXY(0, 0, 1);
    cv::Point3d rotationAxis = normalXY.cross(vector);
    double rotationAngle = RLine::getAngleBetweenVectors(normalXY, vector);
    cv::Mat rotation_1 = RLine::getRotateMatrix(rotationAxis, -rotationAngle);
    cv::Mat rotation_2 = RLine::getRotateMatrix(rotationAxis, rotationAngle);
    cv::Mat rotateVector_1 = rotation_1 * cvPointToMat(vector);
    cv::Mat rotateVector_2 = rotation_2 * cvPointToMat(vector);

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

cv::Point3d RegistrationPrivate::TransformPoint(const cv::Point3d& pPoint, const cv::Mat& pTransform)
{
    cv::Mat pointMat = pTransform * cvPointToMat(pPoint);
    return cv::Point3d(pointMat);
}

double RegistrationPrivate::getACos(double pValue)
{
    if (pValue <= -1.0)
    {
        return acos(-1.0);
    }
    else if (pValue >= 1.0)
    {
        return 0.0;
    }
    else
    {
        return acos(pValue);
    }
}

std::vector<cv::Point3d> RegistrationPrivate::getPointsAtDistanceVTK(vtkIdType centerId, double radius, const vtkSmartPointer<vtkPolyData> pPoly)
{
    pPoly->BuildLinks();
    std::vector<std::pair<cv::Point3d, cv::Point3d>> result;
    std::set<vtkIdType> visited;

    double pnt[3];
    pPoly->GetPoint(centerId, pnt);
    cv::Point3d centerPoint(pnt[0], pnt[1], pnt[2]);

    visited.insert(centerId);
    getPointsAtDistanceRecursiveVTK(centerId, centerPoint, radius, pPoly, result, visited);

    std::vector<cv::Point3d> myPoints;

    for (int i = 0; i < result.size(); i++)
    {
        myPoints.push_back(result[i].second);
    }

    return myPoints;
}

std::pair<double, cv::Point3d> RegistrationPrivate::getMinCircle(const PointTypeITK& a, const PointTypeITK& b, const PointTypeITK& c)
{
    cv::Point3d cvA = itkPointToCV(a);
    cv::Point3d cvB = itkPointToCV(b);
    cv::Point3d cvC = itkPointToCV(c);

    return getMinCircle(cvA, cvB, cvC);
}

std::pair<double, cv::Point3d> RegistrationPrivate::getMinCircle(const cv::Point3d& a, const cv::Point3d& b, const cv::Point3d& c)
{
    std::vector<cv::Point3d> points = { a, b, c };

    RPlane helpPlane;
    helpPlane.init(a, b, c);

    double Z = 0;
    cv::Mat zRot = GetRotateZ(helpPlane.getNormalVector());
    std::vector<cv::Point2f> coplanar2d;

    auto it1 = points.begin();
    auto it2 = points.end();

    for (; it1 != it2; ++it1)
    {
        cv::Point3d pointTemp = TransformPoint(*it1, zRot);
        coplanar2d.push_back(cv::Point2f(pointTemp.x, pointTemp.y));
        Z += pointTemp.z;
    }

    Z = Z / double(points.size());

    cv::Point2f center;
    float radius;

    cv::minEnclosingCircle(coplanar2d, center, radius);

    cv::Point3d tCenter = TransformPoint(cv::Point3d(center.x, center.y, Z), zRot.inv());

    return std::make_pair(radius, tCenter);
}

Eigen::Matrix4d RegistrationPrivate::getRigidTransform(const std::vector<cv::Point3d>& threeVectorsSource, const std::vector<cv::Point3d>& threeVectorstarget, const cv::Point3d& pointSource, const cv::Point3d& pointTarget)
{
    cv::Mat data(3, 1, CV_64F);

    cv::Mat rotation = LeastSquaresICP::GetRotationAnglesXYZ(threeVectorsSource, threeVectorstarget, data);
    cv::Mat translation = cvPointToMat(pointTarget) - (rotation * cvPointToMat(pointSource));

    Eigen::Matrix4d rigidTransform(4, 4);
    rigidTransform(0, 0) = rotation.at<double>(0, 0);
    rigidTransform(1, 0) = rotation.at<double>(1, 0);
    rigidTransform(2, 0) = rotation.at<double>(2, 0);
    rigidTransform(3, 0) = 0;

    rigidTransform(0, 1) = rotation.at<double>(0, 1);
    rigidTransform(1, 1) = rotation.at<double>(1, 1);
    rigidTransform(2, 1) = rotation.at<double>(2, 1);
    rigidTransform(3, 1) = 0;

    rigidTransform(0, 2) = rotation.at<double>(0, 2);
    rigidTransform(1, 2) = rotation.at<double>(1, 2);
    rigidTransform(2, 2) = rotation.at<double>(2, 2);
    rigidTransform(3, 2) = 0;

    rigidTransform(0, 3) = translation.at<double>(0, 0);
    rigidTransform(1, 3) = translation.at<double>(1, 0);
    rigidTransform(2, 3) = translation.at<double>(2, 0);
    rigidTransform(3, 3) = 1;

    return rigidTransform;
}

/*
std::pair<std::vector<cv::Point3d>, std::vector<cv::Point3d>> RegistrationPrivate::getAxisKneeFemurVectors(const PointTypeITK& pHipCamera, const PointTypeITK& pLatEpiCamera, const PointTypeITK& pMedEpiCamera, const PointTypeITK& pHipCT, const PointTypeITK& pLatEpiCT, const PointTypeITK& pMedEpiCT)
{
    /////////////////////////////Camera
    cv::Point3d cameraTEA = itkPointToCV(pLatEpiCamera) - itkPointToCV(pMedEpiCamera);
    RLine cameraLine = RLine::makeLineWithPoints(itkPointToCV(pLatEpiCamera), itkPointToCV(pMedEpiCamera));
    cv::Point3d cameraAxis = itkPointToCV(pHipCamera) - (cameraLine.getProjectPoint(itkPointToCV(pHipCamera)));
    cv::Point3d cameraAP = cameraTEA.cross(cameraAxis);

    std::vector<cv::Point3d> vectorSource = { cameraTEA, cameraAxis, cameraAP };

    ///////////////////////////CT
    cv::Point3d ctTEA = itkPointToCV(pLatEpiCT) - itkPointToCV(pMedEpiCT);
    RLine ctLine = RLine::makeLineWithPoints(itkPointToCV(pLatEpiCT), itkPointToCV(pMedEpiCT));
    cv::Point3d ctAxis = itkPointToCV(pHipCT) - (ctLine.getProjectPoint(itkPointToCV(pHipCT)));
    cv::Point3d ctAP = ctTEA.cross(ctAxis);

    std::vector<cv::Point3d> vectorTarget = { ctTEA, ctAxis, ctAP };

    return std::make_pair(vectorSource, vectorTarget);
}*/

std::pair<std::vector<cv::Point3d>, std::vector<cv::Point3d>> RegistrationPrivate::getAxisPelvisVectors(const cv::Point3d& pAnteriorCamera, const cv::Point3d& pPosteriorCamera, const cv::Point3d& pSuperiorCamera, const cv::Point3d& pAnteriorCT, const cv::Point3d& pPosteriorCT, const cv::Point3d& pSuperiorCT)
{
    ///////////////////////camera

    cv::Point3d posteriorAnteriorCamera = pAnteriorCamera - pPosteriorCamera;
    cv::Point3d posteriorSuperiorCamera = pSuperiorCamera - pPosteriorCamera;
    posteriorAnteriorCamera = posteriorAnteriorCamera / sqrt(posteriorAnteriorCamera.dot(posteriorAnteriorCamera));
    posteriorSuperiorCamera = posteriorSuperiorCamera / sqrt(posteriorSuperiorCamera.dot(posteriorSuperiorCamera));

    cv::Point3d normalCamera = posteriorAnteriorCamera.cross(posteriorSuperiorCamera);
    normalCamera = normalCamera / sqrt(normalCamera.dot(normalCamera));

    posteriorSuperiorCamera = posteriorAnteriorCamera.cross(normalCamera);
    posteriorSuperiorCamera = posteriorSuperiorCamera / sqrt(posteriorSuperiorCamera.dot(posteriorSuperiorCamera));

    std::vector<cv::Point3d> vectorSource = { posteriorAnteriorCamera, posteriorSuperiorCamera, normalCamera };

    //////////////////////////////////CT

    cv::Point3d posteriorAnteriorCT = pAnteriorCT - pPosteriorCT;
    cv::Point3d posteriorSuperiorCT = pSuperiorCT - pPosteriorCT;

    posteriorAnteriorCT = posteriorAnteriorCT / sqrt(posteriorAnteriorCT.dot(posteriorAnteriorCT));
    posteriorSuperiorCT = posteriorSuperiorCT / sqrt(posteriorSuperiorCT.dot(posteriorSuperiorCT));

    cv::Point3d normalCT = posteriorAnteriorCT.cross(posteriorSuperiorCT);
    normalCT = normalCT / sqrt(normalCT.dot(normalCT));

    posteriorSuperiorCT = posteriorAnteriorCT.cross(normalCT);
    posteriorSuperiorCT = posteriorSuperiorCT / sqrt(posteriorSuperiorCT.dot(posteriorSuperiorCT));

    std::vector<cv::Point3d> vectorTarget = { posteriorAnteriorCT, posteriorSuperiorCT, normalCT };

    return std::make_pair(vectorSource, vectorTarget);
}

std::vector<cv::Point3d> RegistrationPrivate::getAxisHipFemoral(const PointTypeITK& pAnteriorFemoralNeck, const PointTypeITK& pAnteriorDistalTrochanter, const PointTypeITK& pLateralTrochanter)
{
    cv::Point3d anteriorFemoralNeckCV, anteriorDistalTrochanterCV, lateralTrochanterCV;
    anteriorFemoralNeckCV = itkPointToCV(pAnteriorFemoralNeck);
    anteriorDistalTrochanterCV = itkPointToCV(pAnteriorDistalTrochanter);
    lateralTrochanterCV = itkPointToCV(pLateralTrochanter);

    return getAxisHipFemoral(anteriorFemoralNeckCV, anteriorDistalTrochanterCV, lateralTrochanterCV);
}

std::vector<cv::Point3d> RegistrationPrivate::getAxisHipFemoral(const cv::Point3d& pAnteriorFemoralNeck, const cv::Point3d& pAnteriorDistalTrochanter, const cv::Point3d& pLateralTrochanter)
{
    cv::Point3d tAxis = pAnteriorFemoralNeck - pAnteriorDistalTrochanter;
    cv::Point3d tVectorHelp = pLateralTrochanter - pAnteriorDistalTrochanter;
    cv::Point3d tNormal = tAxis.cross(tVectorHelp);
    cv::Point3d tCross = tAxis.cross(tNormal);

    cv::Point3d VectorAxisTop = tAxis / sqrt(tAxis.dot(tAxis));
    cv::Point3d VectorNormal = tNormal / sqrt(tNormal.dot(tNormal));
    cv::Point3d CrossVector = tCross / sqrt(tCross.dot(tCross));

    return { VectorAxisTop, VectorNormal, CrossVector };
}

cv::Mat RegistrationPrivate::GetRotateTransformVectors(const cv::Point3d& pFromVector, const cv::Point3d& pToVector)
{
    cv::Point3d rotationAxis = pFromVector.cross(pToVector);

    double rotationAngle = getAngleBetweenVectors(pFromVector, pToVector);

    cv::Mat rotate = getRotateMatrix(rotationAxis, rotationAngle);

    return rotate;
}

double RegistrationPrivate::getAngleBetweenVectors(const cv::Point3d& a, const cv::Point3d& b)
{
    double scalar = a.dot(b);
    double magnitude = sqrt((a.dot(a)) * (b.dot(b)));
    double tCos = scalar / magnitude;
    if (tCos <= -1.0)
    {
        return EIGEN_PI;
    }
    else if (tCos >= 1.0)
    {
        return 0;
    }
    else
    {
        return acos(tCos);
    }
}

cv::Mat RegistrationPrivate::getRotateMatrix(const cv::Point3d& axis, double angle)
{
    cv::Mat rotationMatrix(3, 3, CV_64F);

    cv::Point3d normaliceAxis = axis;
    
    normaliceAxis = normaliceAxis / sqrt(normaliceAxis.dot(normaliceAxis));

    rotationMatrix.at <double>(0, 0) = 0;
    rotationMatrix.at <double>(1, 0) = normaliceAxis.z;
    rotationMatrix.at <double>(2, 0) = -normaliceAxis.y;

    rotationMatrix.at <double>(0, 1) = -normaliceAxis.z;
    rotationMatrix.at <double>(1, 1) = 0;
    rotationMatrix.at <double>(2, 1) = normaliceAxis.x;

    rotationMatrix.at <double>(0, 2) = normaliceAxis.y;
    rotationMatrix.at <double>(1, 2) = -normaliceAxis.x;
    rotationMatrix.at <double>(2, 2) = 0;

    cv::Mat I = cv::Mat::eye(3, 3, CV_64F);
    cv::Mat result = I + sin(angle)*rotationMatrix + (1.0 - cos(angle))*(rotationMatrix * rotationMatrix);
    return result;
}

void RegistrationPrivate::getPointsAtDistanceRecursiveVTK(vtkIdType baseId, const cv::Point3d& pCenterPoint, double radius, const vtkSmartPointer<vtkPolyData> pPoly, std::vector<std::pair<cv::Point3d, cv::Point3d>>& result, std::set<vtkIdType>& visited)
{
    if (result.size() > 1)
    {
        return;
    }

    double pnt[3];
    pPoly->GetPoint(baseId, pnt);
    cv::Point3d basePoint(pnt[0], pnt[1], pnt[2]);

    vtkNew<vtkIdList> cells;
    pPoly->GetPointCells(baseId, cells);

    vtkIdType nCells = cells->GetNumberOfIds();

    for (int i = 0; i < nCells; i++)
    {
        vtkIdType cellId = cells->GetId(i);

        vtkSmartPointer<vtkCell> cell = pPoly->GetCell(cellId);
        vtkSmartPointer<vtkIdList> cellPoints = cell->GetPointIds();

        for (int j = 0; j < cellPoints->GetNumberOfIds(); j++)
        {
            vtkIdType pointTempId = cellPoints->GetId(j);
            if (visited.find(pointTempId) == visited.end())
            {
                visited.insert(pointTempId);

                double pnt[3];
                pPoly->GetPoint(pointTempId, pnt);
                cv::Point3d pointTemp(pnt[0], pnt[1], pnt[2]);
                
                if (GetDistance(pointTemp, pCenterPoint) >= radius)
                {
                    result.push_back(std::make_pair(basePoint, pointTemp));
                }
                else
                {
                    getPointsAtDistanceRecursiveVTK(pointTempId, pCenterPoint, radius, pPoly, result, visited);
                }
            }
        }
    }
}
