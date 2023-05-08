#include "Segment.hpp"
#include "Utils.hpp"

using namespace THA::IMPLANTS;

JoinedSegment::JoinedSegment(const SegmentPoint2D& pA, const SegmentPoint2D& pB, const SegmentPoint2D& pC)
{
    A_2D = pA;
    B_2D = pB;
    C_2D = pC;
    is2D = true;
}

JoinedSegment::JoinedSegment(const SegmentPoint3D& pA, const SegmentPoint3D& pB, const SegmentPoint3D& pC)
{
    A_3D = pA;
    B_3D = pB;
    C_3D = pC;
    is2D = false;
}

double JoinedSegment::getAngle() const
{
    double scalar, magnitude;

    if (is2D == true)
    {
        cv::Point2f a = A_2D.point - B_2D.point;
        cv::Point2f b = C_2D.point - B_2D.point;
        scalar = a.dot(b);
        magnitude = sqrt((a.dot(a)) * (b.dot(b)));
    }
    else
    {
        Point a = A_3D.point - B_3D.point;
        Point b = C_3D.point - B_3D.point;
        scalar = a.dot(b);
        magnitude = sqrt((a.dot(a)) * (b.dot(b)));
    }

    double tCos = scalar / magnitude;
    if (tCos <= -1.0)
    {
        return PI;
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

bool JoinedSegment::operator< (const JoinedSegment& segment) const
{
    return (this->getAngle() > segment.getAngle());
}

int JoinedSegment::GetCenterPos() const
{
    if (is2D == true)
    {
        return B_2D.pos;
    }
    else
    {
        return B_3D.pos;
    }
}

//cv::Point2f JoinedSegment2D::GetA() const
//{
//    return A;
//}
//
//cv::Point2f JoinedSegment2D::GetC() const
//{
//    return C;
//}
//
//cv::Point2f JoinedSegment2D::GetO() const
//{
//    return B;
//}
//
//float JoinedSegment2D::GetDistance() const
//{
//    cv::Point2f diff = A - C;
//    float square = diff.dot(diff);
//    return sqrt(square);
//}
//
//cv::Point2f JoinedSegment2D::GetCenterPoint() const
//{
//    return (A + C) / 2.0;
//}