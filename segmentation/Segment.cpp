#include "Segment.hpp"


JoinedSegment::JoinedSegment(const SegmentPoint& pA, const SegmentPoint& pB, const SegmentPoint& pC)
{
    A_2D = pA;
    B_2D = pB;
    C_2D = pC;
}

double JoinedSegment::getAngle() const
{
    double scalar, magnitude;

    cv::Point3f a = A_2D.point - B_2D.point;
    cv::Point3f b = C_2D.point - B_2D.point;
    scalar = a.dot(b);
    magnitude = sqrt((a.dot(a)) * (b.dot(b)));

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
   return B_2D.pos;
}
