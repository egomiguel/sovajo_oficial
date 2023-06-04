#ifndef THA_SEGMENT2D_IMPLANT_H
#define THA_SEGMENT2D_IMPLANT_H

#include <opencv2/calib3d/calib3d.hpp>
#include "Point.hpp"

namespace THA
{
	namespace IMPLANTS
	{
		struct SegmentPoint2D
		{
			cv::Point2f point;
			int pos;
		};

		struct SegmentPoint3D
		{
			Point point;
			int pos;
		};

		class JoinedSegment
		{
		private:
			SegmentPoint2D A_2D, B_2D, C_2D;
			SegmentPoint3D A_3D, B_3D, C_3D;
			bool is2D;
		public:
			JoinedSegment(const SegmentPoint2D& pA, const SegmentPoint2D& pB, const SegmentPoint2D& pC);
			JoinedSegment(const SegmentPoint3D& pA, const SegmentPoint3D& pB, const SegmentPoint3D& pC);
			double getAngle() const;
			bool operator< (const JoinedSegment& segment) const;
			int GetCenterPos() const;
			/*cv::Point2f GetA() const;
			cv::Point2f GetC() const;
			cv::Point2f GetO() const;
			cv::Point2f GetCenterPoint() const;
			float GetDistance() const;*/
		};
	}
}

#endif