#include <iostream>
#include <math.h>
#include <numeric>
#include "DataFit.hpp"

using namespace TKA::HIP;

DataFit::DataFit(const std::vector<cv::Point3d>& ellipse_points)
{
	ellipse_points_ = ellipse_points;
	centroid_ = GetCentroid(ellipse_points_);
	plane_ = FitPlane(ellipse_points_, centroid_);

	//b*X - a * Y + 0 = b * X1 - a * Y1
	normal_.equation1.x = plane_[1];
	normal_.equation1.y = -plane_[0];
	normal_.equation1.z = 0.0;
	normal_.bias1 = plane_[1] * centroid_.x - plane_[0] * centroid_.y;

	//c*X + 0 - a * Z = c * X1 - a * Z1
	normal_.equation2.x = plane_[2];
	normal_.equation2.y = 0.0;
	normal_.equation2.z = -plane_[0];
	normal_.bias2 = plane_[2] * centroid_.x - plane_[0] * centroid_.z;
}

cv::Point3d DataFit::GetCentroid(const std::vector<cv::Point3d>& points)
{
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;

	for (int i = 0; i < points.size(); i++)
	{
		if (i == 0)
		{
			x.push_back(points[i].x);
			x.push_back(points[i].x);
			x.push_back(points[i].x);
			x.push_back(points[i].x);

			y.push_back(points[i].y);
			y.push_back(points[i].y);
			y.push_back(points[i].y);
			y.push_back(points[i].y);

			z.push_back(points[i].z);
			z.push_back(points[i].z);
			z.push_back(points[i].z);
			z.push_back(points[i].z);
		}

		if (points[i].x > x[0])
		{
			x[1] = x[0];
			x[0] = points[i].x;
		}
		else
		{
			if (points[i].x < x[2])
			{
				x[3] = x[2];
				x[2] = points[i].x;

			}
		}

		if (points[i].y > y[0])
		{
			y[1] = y[0];
			y[0] = points[i].y;
		}
		else
		{
			if (points[i].y < y[2])
			{
				y[3] = y[2];
				y[2] = points[i].y;

			}
		}

		if (points[i].z > z[0])
		{
			z[1] = z[0];
			z[0] = points[i].z;
		}
		else
		{
			if (points[i].z < z[2])
			{
				z[3] = z[2];
				z[2] = points[i].z;

			}
		}
	}

	return cv::Point3d(std::accumulate(x.begin(), x.end(), 0) / x.size(), std::accumulate(y.begin(), y.end(), 0) / y.size(), std::accumulate(z.begin(), z.end(), 0) / z.size());
}

std::vector<double> DataFit::FitPlane(std::vector<cv::Point3d>& points, const cv::Point3d& centroid)
{
	//a*x + b*y + c*z = d
	int rows = points.size();
	cv::Mat A = cv::Mat(rows, 3, CV_64F, points.data());
	//Subtracting centroid_ from the set of points to center on the origin.

    std::vector<cv::Point3d> center(rows, centroid);

    cv::Mat C(rows, 3, CV_64F, center.data());

	cv::Mat CA = A - C;
	cv::SVD svd(CA);
	cv::Mat Last_V_Row = svd.vt.row(2);
	double a = Last_V_Row.at<double>(0);
	double b = Last_V_Row.at<double>(1);
	double c = Last_V_Row.at<double>(2);
	double d = a * centroid.x + b * centroid.y + c * centroid.z;
	std::vector<double> result{ a, b, c, d };
	return result;
}

Normal DataFit::GetNormalEquation()
{
	return normal_;
}