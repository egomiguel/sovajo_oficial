#include "SliceBorderPrivate.hpp"



SliceBorderPrivate::SliceBorderPrivate()
{
}

double SliceBorderPrivate::GetDistanceBetweenPoints(cv::Point3d a, cv::Point3d b)
{
    cv::Point3d diff = a - b;
    return diff.dot(diff);
}

bool SliceBorderPrivate::DeletePointsInsideRadius(std::list<cv::Point3d>& points, const cv::Point3d& centerPoint, cv::Point3d& nearPoint, double radius)
{
    bool result = false;
    double distance;
    double closeDistance = 9999999.0;
    std::list<cv::Point3d>::iterator nearPointIt;
    auto it = points.begin();
    while (it != points.end())
    {
        distance = GetDistanceBetweenPoints(centerPoint, *it);
        if (distance < radius)
        {
            it = points.erase(it);
        }
        else
        {
            if (distance < closeDistance)
            {
                closeDistance = distance;
                nearPointIt = it;
                nearPoint = *it;
                result = true;
            }
            ++it;
        }
    }
    if (result == true)
    {
        points.erase(nearPointIt);
    }
    return result;
}

//std::vector<cv::Point3d> SliceBorderPrivate::SortPoints(const std::list<cv::Point3d>& points)
//{
//    cv::Point3d nearPoint, changePoint;
//    std::vector<cv::Point3d> result;
//    std::list<cv::Point3d> resultTemp = points;
//    if (resultTemp.size() > 0)
//    {
//        changePoint = *(resultTemp.begin());
//        resultTemp.erase(resultTemp.begin());
//        result.push_back(changePoint);
//    }
//    else
//    {
//        return result;
//    }
//
//    std::list<cv::Point3d>::iterator itList1, itTemp;
//
//    do
//    {
//        if (DeletePointsInsideRadius(resultTemp, changePoint, nearPoint, 0.2))
//        {
//            changePoint = nearPoint;
//            result.push_back(nearPoint);
//        }
//
//    } while (resultTemp.size() > 0);
//
//    return result;
//}

int SliceBorderPrivate::GetLessAnglePosition(const std::vector<JoinedSegment>& segments, double& angle)
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
    angle = (*lessAngle).getAngle();
    return (*lessAngle).GetCenterPos();
}

std::vector<cv::Point3d> SliceBorderPrivate::GetContourMainPoints(const std::vector<cv::Point3d>& sortPoints, double angle)
{
    std::vector<JoinedSegment> joinedSegment;
    std::list<std::pair<cv::Point3d, int> > data;
    std::list<std::pair<cv::Point3d, int> >::iterator moveIt, it1, it2, initIt, endIt;

    std::vector<cv::Point3d> result = sortPoints;

    for (int i = 0; i < result.size(); i++)
    {
        data.push_back(std::make_pair(result[i], i));
    }

    if (result.size() == 0)
    {
        return result;
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
                SegmentPoint a, b, c;

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
                SegmentPoint a, b, c;

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
                SegmentPoint a, b, c;

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
        double myAngle;
        double radAngle = angle * PI / 180.0;

        int pos = GetLessAnglePosition(joinedSegment, myAngle);
        if (myAngle >= radAngle)
        {
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
        }
        else
        {
            break;
        }

    } while (data.size() > 0);

    it1 = data.begin();
    it2 = data.end();

    for (; it1 != it2; ++it1)
    {
        result.push_back((*it1).first);
    }

    return result;

}

std::vector<cv::Point3d> SliceBorderPrivate::InterpolateSpline(const std::vector<cv::Point3d>& pPoints, int amount)
{
    std::vector<cv::Point3d> result;

    if (pPoints.size() == 0)
    {
        return result;
    }

    Eigen::MatrixXd points(3, pPoints.size());
    int row_index = 0;
    for (auto const point : pPoints)
    {
        points.col(row_index) << point.x, point.y, point.z;
        row_index++;
    }
    Eigen::Spline3d spline = Eigen::SplineFitting<Eigen::Spline3d>::Interpolate(points, 2);

    double pos = 0;
    for (int i = 0; i < amount; i++) {
        pos += 1.0 / (amount * 1.0);
        Eigen::VectorXd values = spline(pos);
        result.push_back(cv::Point3d(values[0], values[1], values[2]));
    }

    return result;
}

cv::Point3d SliceBorderPrivate::ArrayToPoint(const double point[3])
{
    cv::Point3d result = { point[0], point[1], point[2] };
    return result;
}