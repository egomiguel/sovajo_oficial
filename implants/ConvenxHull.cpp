#include "ConvexHull.hpp"
#include <list>
#include <iostream>
#include "opencv2/imgproc.hpp"
//#include <pcl/surface/concave_hull.h>
#include <Eigen/Core>
#include <unsupported/Eigen/Splines>
#include "Line.hpp"

struct GeneralLine
{
    float A, B, C;
    GeneralLine(const cv::Point2f& p1, const cv::Point2f& p2)
    {
		A = -(p2.y - p1.y);
		B = p2.x - p1.x;
		C = -(A * p1.x + B * p1.y);

		float norm = sqrt(A * A + B * B);
		A = A / norm;
		B = B / norm;
		C = C / norm;
    }

    float eval(const cv::Point2f& p)
    {
        return (A * p.x + B * p.y + C);
    }

};

ConvexHull::ConvexHull(const std::vector<Point>& pPoints, const cv::Mat& pRotationOnZ)
{
    rotationOnZInv = pRotationOnZ.inv();
	rotationOnZ = pRotationOnZ;
    axisZ = 0;
    auto it1 = pPoints.begin();
    auto it2 = pPoints.end();

    for (; it1 != it2; ++it1)
    {
        cv::Mat rotatePointMat = rotationOnZ * (*it1).ToMatPoint();
        cv::Point3d rotatePoint = cv::Point3d(rotatePointMat);
        cv::Point3f fPoint = static_cast<cv::Point3f>(rotatePoint);
        points.push_back(cv::Point2f(fPoint.x, fPoint.y));
        axisZ += fPoint.z;
    }

    if (pPoints.size() > 0)
    {
        axisZ = axisZ / float(pPoints.size());
    }

	cv::convexHull(points, mConvexHull2D, true);

	auto hullIt1 = mConvexHull2D.begin();
	auto hullIt2 = mConvexHull2D.end();

	mCenterPoint = cv::Point2f(0, 0);

	for (; hullIt1 != hullIt2; ++hullIt1)
	{
		mCenterPoint = mCenterPoint + (*hullIt1);
		cv::Point2d tempDouble = static_cast<cv::Point2d>(*hullIt1);
		Point tempPoint(tempDouble.x, tempDouble.y, static_cast<double>(axisZ));
		cv::Mat rotatePointMat = rotationOnZInv * tempPoint.ToMatPoint();
		mConvexHull.push_back(Point(rotatePointMat));
	}

	if (mConvexHull2D.size() > 0)
	{
		mCenterPoint = mCenterPoint / float(mConvexHull2D.size());
	}

}

bool ConvexHull::isPointWithinConvexHull(const Point& pPoint, float pErrorMargin)
{
	/*
		The points are supposed to be clockwise, therefore a point that is inside 
		the convex hull will always be to the right of each segment and its 
		evaluation will be negative.
	*/

	cv::Mat rotatePointMat = rotationOnZ * pPoint.ToMatPoint();
	Point rotationPoint = Point(rotatePointMat);
	std::vector<cv::Point2f> tempConvex;

	if (pErrorMargin <= 0)
	{
		tempConvex = mConvexHull2D;
	}
	else
	{
		tempConvex = increaseConvexHull(pErrorMargin);
	}

	int tSize = tempConvex.size();
	cv::Point2d rotationPoint2D(rotationPoint.x, rotationPoint.y);

	for (int i = 0; i < tSize; i++)
	{
		cv::Point2d p1 = tempConvex[i];
		cv::Point2d p2 = tempConvex[(i + 1) % tSize];
		GeneralLine tempLine(p1, p2);
		float result = tempLine.eval(rotationPoint2D);
		//std::cout << pErrorMargin << " Dista: " << result << std::endl;
		if (result > 0)
		{
			return false;
		}
	}
	return true;
}

bool ConvexHull::areSomePointWithinConvexHull(const std::vector<Point>& pPoints, float pErrorMargin)
{
	int tSize = pPoints.size();
	for (int i = 0; i < tSize; i++)
	{
		bool result = isPointWithinConvexHull(pPoints[i], pErrorMargin);

		if (result == true)
		{
			return true;
		}
	}
	return false;
}

std::vector<cv::Point2f> ConvexHull::increaseConvexHull(float increaseDist)
{
	int tSize = mConvexHull2D.size();
	std::vector<cv::Point2f> increasePoints, newConvexHull;

	for (int i = 0; i < tSize; i++)
	{
		cv::Point2f vector = mConvexHull2D[i] - mCenterPoint;
		vector = vector / sqrt(vector.dot(vector));
		increasePoints.push_back(mConvexHull2D[i] + increaseDist * vector);
	}

	cv::convexHull(increasePoints, newConvexHull, true);
	return newConvexHull;
}

std::vector<Point> ConvexHull::GetConvexHull(int vertices)
{

    if (vertices < 3 || vertices >= mConvexHull.size())
    {
        return mConvexHull;
    }

	std::vector<Point> result = mConvexHull;

    std::vector<JoinedSegment> joinedSegment;
    std::list<std::pair<Point, int> > data;
    std::list<std::pair<Point, int> >::iterator moveIt, it1, it2, initIt, endIt;

    for (int i = 0; i < result.size(); i++)
    {
        data.push_back(std::make_pair(result[i], i));
    }
    result.clear();
    do
    {
        initIt = data.begin();
        endIt = data.end();
        it1 = std::next(initIt, data.size() - 1);
        it2 = std::next(initIt, data.size() - 2);

        for (; initIt != endIt; ++initIt)
        {
            if (initIt == it2)
            {
                SegmentPoint3D a, b, c;

                moveIt = initIt;

                a.point = (*moveIt).first;
                a.pos = (*moveIt).second;

                moveIt = std::next(initIt, 1);

                b.point = (*moveIt).first;
                b.pos = (*moveIt).second;

                moveIt = data.begin();

                c.point = (*moveIt).first;
                c.pos = (*moveIt).second;

                joinedSegment.push_back(JoinedSegment(a, b, c));
            }
            else if (initIt == it1)
            {
                SegmentPoint3D a, b, c;

                moveIt = initIt;

                a.point = (*moveIt).first;
                a.pos = (*moveIt).second;

                moveIt = data.begin();

                b.point = (*moveIt).first;
                b.pos = (*moveIt).second;

                moveIt = std::next(data.begin(), 1);

                c.point = (*moveIt).first;
                c.pos = (*moveIt).second;

                joinedSegment.push_back(JoinedSegment(a, b, c));
            }
            else
            {
                SegmentPoint3D a, b, c;

                moveIt = initIt;

                a.point = (*moveIt).first;
                a.pos = (*moveIt).second;

                moveIt = std::next(initIt, 1);

                b.point = (*moveIt).first;
                b.pos = (*moveIt).second;

                moveIt = std::next(initIt, 2);

                c.point = (*moveIt).first;
                c.pos = (*moveIt).second;

                joinedSegment.push_back(JoinedSegment(a, b, c));
            }
        }

        int pos = GetLessAnglePosition(joinedSegment);
        joinedSegment.clear();

        it1 = data.begin();
        it2 = data.end();

        for (; it1 != it2; ++it1)
        {
            if ((*it1).second == pos)
            {
                data.erase(it1);
                break;
            }
        }

    } while (data.size() > vertices);

    it1 = data.begin();
    it2 = data.end();

    for (; it1 != it2; ++it1)
    {
        result.push_back((*it1).first);
    }
    
    return result;

}

std::vector<Point> ConvexHull::ReduceConvexHull(const std::vector<Point>& pPoints, int vertices)
{

    if (vertices < 3 || vertices >= pPoints.size())
    {
        return pPoints;
    }

    std::vector<JoinedSegment> joinedSegment;
    std::list<std::pair<Point, int> > data;
    std::list<std::pair<Point, int> >::iterator moveIt, it1, it2;

    for (int i = 0; i < pPoints.size(); i++)
    {
        data.push_back(std::make_pair(pPoints[i], i));
    }

    do
    {
        it1 = data.begin();
        it2 = std::next(it1, data.size() - 2);

        for (; it1 != it2; ++it1)
        {
            SegmentPoint3D a, b, c;

            moveIt = it1;

            a.point = (*moveIt).first;
            a.pos = (*moveIt).second;

            moveIt = std::next(it1, 1);

            b.point = (*moveIt).first;
            b.pos = (*moveIt).second;

            moveIt = std::next(it1, 1);

            c.point = (*moveIt).first;
            c.pos = (*moveIt).second;

            joinedSegment.push_back(JoinedSegment(a, b, c));
        }

        int pos = GetLessAnglePosition(joinedSegment);
        joinedSegment.clear();

        it1 = data.begin();
        it2 = data.end();

        for (; it1 != it2; ++it1)
        {
            if ((*it1).second == pos)
            {
                data.erase(it1);
                break;
            }
        }

    } while (data.size() > vertices);

    std::vector<Point> result;

    it1 = data.begin();
    it2 = data.end();

    for (; it1 != it2; ++it1)
    {
        result.push_back((*it1).first);
    }

    return result;
}

std::vector<Point> ConvexHull::interpolateSpline(const std::vector<Point>& pPoints, int amount)
{
    std::vector<Point> result;

    Eigen::MatrixXd points(3, pPoints.size());
    int row_index = 0;
    for (auto const way_point : pPoints)
    {
        points.col(row_index) << way_point.x, way_point.y, way_point.z;
        row_index++;
    }
    Eigen::Spline3d spline = Eigen::SplineFitting<Eigen::Spline3d>::Interpolate(points, 2);

    double time_ = 0;
    for (int i = 0; i < amount; i++) {
        time_ += 1.0 / (amount * 1.0);
        Eigen::VectorXd values = spline(time_);
        result.push_back(Point(values[0], values[1], values[2]));
    }

    return result;
}

int ConvexHull::GetLessAnglePosition(const std::vector<JoinedSegment>& segments)
{
    std::vector<JoinedSegment>::const_iterator lessAngle = segments.begin();

    auto it1 = segments.begin();
    auto it2 = segments.end();

    for (; it1 != it2; ++it1)
    {
        if ((*it1) < (*lessAngle))
        {
            lessAngle = it1;
        }
    }
    return (*lessAngle).GetCenterPos();
}

//
//ConvexHull::ConvexHull(const std::vector<cv::Point2d>& pPoints)
//{
//    auto it1 = pPoints.begin();
//    auto it2 = pPoints.end();
//    for (; it1 != it2; ++it1)
//    {
//        points.push_back(static_cast<cv::Point2f>(*it1));
//    }
//}
//
//std::vector<cv::Point2d> ConvexHull::GetPointsDouble(const std::vector<cv::Point2f>& pPoints)
//{
//    std::vector<cv::Point2d> result;
//    auto it1 = pPoints.begin();
//    auto it2 = pPoints.end();
//    for (; it1 != it2; ++it1)
//    {
//        result.push_back(static_cast<cv::Point2d>(*it1));
//    }
//    return result;
//}
//
//// 3D cross product of OA and OB vectors, (i.e z - component of their "2D" cross product, but remember that it is not defined in "2D").
//// Returns a positive value, if OAB makes a counter-clockwise turn,
//// negative for clockwise turn, and zero if the points are collinear.
//int ConvexHull::orientation(const cv::Point2d& O, const cv::Point2d& A, const cv::Point2d& B)
//{
//    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
//}
//
//std::vector<cv::Point2d> ConvexHull::GetConvexHull(int vertices)
//{
//    std::vector<cv::Point2f> result;
//    cv::convexHull(points, result, true);
//
//    if (vertices < 3 || vertices >= result.size())
//    {
//        return GetPointsDouble(result);
//    }
//
//    std::vector<JoinedSegment> joinedSegment;
//    std::list<std::pair<cv::Point2f, int> > data;
//    std::list<std::pair<cv::Point2f, int> >::iterator moveIt, it1, it2;
//
//    for (int i = 0; i < result.size(); i++)
//    {
//        data.push_back(std::make_pair(result[i], i));
//    }
//
//    do
//    {
//        it1 = data.begin();
//        it2 = std::next(it1, data.size() - 2);
//
//        for ( ; it1 != it2; ++it1 )
//        {
//            /*
//            if (i == data.size() - 2)
//            {
//                SegmentPoint2D a, b, c;
//                a.point = data[i].first;
//                a.pos = data[i].second;
//
//                b.point = data[i + 1].first;
//                b.pos = data[i + 1].second;
//
//                c.point = data[0].first;
//                c.pos = data[0].second;
//
//                joinedSegment.push_back(JoinedSegment(a, b, c));
//            }
//            else if (i == data.size() - 1)
//            {
//                SegmentPoint2D a, b, c;
//
//                a.point = data[i].first;
//                a.pos = data[i].second;
//
//                b.point = data[0].first;
//                b.pos = data[0].second;
//
//                c.point = data[1].first;
//                c.pos = data[1].second;
//
//                joinedSegment.push_back(JoinedSegment(a, b, c));
//            }*/
//
//            SegmentPoint2D a, b, c;
//
//            moveIt = it1;
//
//            a.point = (*moveIt).first;
//            a.pos = (*moveIt).second;
//
//            moveIt = std::next(it1, 1);
//
//            b.point = (*moveIt).first;
//            b.pos = (*moveIt).second;
//
//            moveIt = std::next(it1, 1);
//
//            c.point = (*moveIt).first;
//            c.pos = (*moveIt).second;
//
//            joinedSegment.push_back(JoinedSegment(a, b, c));
//        }
//
//        int pos = GetLessAnglePosition(joinedSegment);
//        joinedSegment.clear();
//
//        it1 = data.begin();
//        it2 = data.end();
//
//        for (; it1 != it2; ++it1)
//        {
//            if ((*it1).second == pos)
//            {
//                data.erase(it1);
//                break;
//            }
//        }
//
//    } while (data.size() > vertices);
//
//    std::vector<cv::Point2f> temp;
//
//    it1 = data.begin();
//    it2 = data.end();
//
//    for ( ; it1 != it2; ++it1 )
//    {
//        temp.push_back((*it1).first);
//    }
//    
//    return GetPointsDouble(temp);
//}

/*
std::vector<Point> ConvexHull::GetConcaveHull(std::vector<Point> coplanar)
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr coplanarPCL(new pcl::PointCloud<pcl::PointXYZ>);

    auto it1 = coplanar.begin();
    auto it2 = coplanar.end();

    for (; it1 != it2; ++it1)
    {
        coplanarPCL->points.push_back(pcl::PointXYZ((*it1).x, (*it1).y, (*it1).z));
    }

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_hull(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::ConcaveHull<pcl::PointXYZ> chull;
    chull.setInputCloud(coplanarPCL);
    chull.setAlpha(0.5);
    chull.reconstruct(*cloud_hull);
    std::vector<Point> result;

    auto pclIt1 = cloud_hull->points.begin();
    auto pclIt2 = cloud_hull->points.end();

    for (; pclIt1 != pclIt2; ++pclIt1)
    {
        result.push_back(Point((*pclIt1).x, (*pclIt1).y, (*pclIt1).z));
    }

    return result;
}*/


