#ifndef SEGMENT2D_IMPLANT_H
#define SEGMENT2D_IMPLANT_H

#include <opencv2/calib3d/calib3d.hpp>

namespace TKA
{
	namespace SEGMENTATION
	{

		const double PI = acos(-1.0);

		struct SegmentPoint
		{
			cv::Point3f point;
			int pos;
		};

		class JoinedSegment
		{
		private:
			SegmentPoint A_2D, B_2D, C_2D;
		public:
			JoinedSegment(const SegmentPoint& pA, const SegmentPoint& pB, const SegmentPoint& pC);
			double getAngle() const;
			bool operator< (const JoinedSegment& segment) const;
			int GetCenterPos() const;
		};
	}
}

#endif