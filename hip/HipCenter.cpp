#include <iostream>
#include <math.h>
#include <numeric>
#include <fstream>
#include "HipCenter.hpp"
#include "HipException.hpp"

using namespace TKA::HIP;

HipCenter::HipCenter(const std::vector<std::vector<HipPoint> >& ellipse_list)
{
	is_ellipse_ = true;
	int tSize = ellipse_list.size();
	if ( tSize < 3)
	{
        throw HipExceptionCode::YOU_HAVE_NOT_GENERATED_ENOUGH_DATA;
	}

	for (int i = 0; i < tSize; i++)
	{
		for (int j = 0; j < ellipse_list[i].size(); j++)
		{
			cv::Point3d tPoint(ellipse_list[i][j].x, ellipse_list[i][j].y, ellipse_list[i][j].z);
			ellipse_list_[i].push_back(tPoint);
		}
	}
}

HipCenter::HipCenter(const std::vector<HipPoint>& point_list)
{
	is_ellipse_ = false;
	int tSize = point_list.size();
	if (tSize < 20)
	{
        throw HipExceptionCode::YOU_HAVE_NOT_GENERATED_ENOUGH_DATA;
	}

	for (int i = 0; i < tSize; i++)
	{
		cv::Point3d tPoint(point_list[i].x, point_list[i].y, point_list[i].z);
		point_list_.push_back(tPoint);
	}
}

cv::Point3d HipCenter::LeastSquareSolve(const cv::Mat& A, const cv::Mat& B) const
{
	cv::SVD svd(A);
	cv::Mat pinvA = svd.vt.t()* cv::Mat::diag(1. / svd.w)*svd.u.t();
	cv::Mat X = pinvA * cv::Mat(B);
	return cv::Point3d(X.at<double>(0, 0), X.at<double>(1, 0), X.at<double>(2, 0));
}

void HipCenter::GetHipCenterByEllipses(double point_out[3]) const
{
	if (is_ellipse_ == false)
	{
        throw HipExceptionCode::ELLIPSES_HAVE_NOT_BEEN_INITIALIZED;
	}

	std::vector<std::vector<cv::Point3d>>::const_iterator it1, it2;
	it1 = ellipse_list_.begin();
	it2 = ellipse_list_.end();
	Normal normal;
	std::vector<cv::Point3d> equations;
	std::vector<double> bias;
	for (; it1 != it2; ++it1)
	{
		DataFit fit(*it1);
		normal = fit.GetNormalEquation();
		equations.push_back(normal.equation1);
		equations.push_back(normal.equation2);
		bias.push_back(normal.bias1);
		bias.push_back(normal.bias2);
	}

	cv::Mat A = cv::Mat(equations.size(), 3, CV_64F, equations.data());
	cv::Mat B = cv::Mat(bias.size(), 1, CV_64F, bias.data());
	cv::Point3d X = LeastSquareSolve(A, B);

	point_out[0] = X.x;
	point_out[1] = X.y;
	point_out[2] = X.z;
}

void HipCenter::GetHipCenterBySphere(double point_out[3]) const
{
	if (is_ellipse_ == true)
	{
        throw HipExceptionCode::SPHERE_HAVE_NOT_BEEN_INITIALIZED;
	}

	cv::Point3d X = GetSphereCenter(point_list_);

	point_out[0] = X.x;
	point_out[1] = X.y;
	point_out[2] = X.z;
}

void HipCenter::GetHipCenterByThreeSphere(double point_out[3]) const
{
	int tSize = point_list_.size();
	if (is_ellipse_ == true)
	{
        throw HipExceptionCode::SPHERE_HAVE_NOT_BEEN_INITIALIZED;
	}

	if (tSize < 60)
	{
        throw HipExceptionCode::YOU_HAVE_NOT_GENERATED_ENOUGH_DATA;
	}

	std::vector<cv::Point3d> sphere_1, sphere_2, sphere_3;
	
	for (int i = 0; i < tSize; i++)
	{
		if (i % 3 == 0)
		{
			sphere_1.push_back(point_list_[i]);
		}
		else if ((i - 1) % 3 == 0)
		{
			sphere_2.push_back(point_list_[i]);
		}
		else
		{
			sphere_3.push_back(point_list_[i]);
		}
	}

	cv::Point3d X = GetSphereCenter(sphere_1);
	cv::Point3d Y = GetSphereCenter(sphere_2);
	cv::Point3d Z = GetSphereCenter(sphere_3);

	cv::Point3d diff1 = X - Y;
	cv::Point3d diff2 = X - Z;
	cv::Point3d diff3 = Y - Z;

	double distance_xy = diff1.dot(diff1);
	double distance_xz = diff2.dot(diff2);
	double distance_yz = diff3.dot(diff3);

	std::vector<double> distances = { distance_xy , distance_xz , distance_yz };
	std::sort(distances.begin(), distances.end());
	
	cv::Point3d result;

	if (distances[0] == distance_xy)
	{
		result = (X + Y) / 2.0;
	}
	else if (distances[0] == distance_xz)
	{
		result = (X + Z) / 2.0;
	}
	else
	{
		result = (Y + Z) / 2.0;
	}

	point_out[0] = result.x;
	point_out[1] = result.y;
	point_out[2] = result.z;
}

cv::Point3d HipCenter::GetSphereCenter(const std::vector<cv::Point3d>& sphere_points) const
{
	std::vector<cv::Point3d>::const_iterator it1, it2;
	it1 = sphere_points.begin();
	it2 = sphere_points.end();
	std::vector<cv::Mat> points;
	std::vector<double> bias;

	cv::Mat A(sphere_points.size(), 4, CV_64F);
	int cont = 0;
	for (; it1 != it2; ++it1)
	{
		cv::Point3d temp = (*it1);
		A.at<double>(cont, 0) = temp.x;
		A.at<double>(cont, 1) = temp.y;
		A.at<double>(cont, 2) = temp.z;
		A.at<double>(cont, 3) = 1.0;
		bias.push_back(temp.dot(temp));
		cont++;
	}

	cv::Mat B = cv::Mat(bias.size(), 1, CV_64F, bias.data());
	cv::Point3d X = LeastSquareSolve(A, B);

	X = X / 2.0;

	return X;
}

void HipCenter::GetHipCenterFromFile(std::string path, double point_out[3])
{
	std::ifstream infile(path);
	std::vector<HipPoint> points;
	double a, b, c;
	int cont = 0;
	HipPoint myPoint;
	while (infile >> a >> b >> c)
	{
		myPoint.x = a;
		myPoint.y = b;
		myPoint.z = c;
		points.push_back(myPoint);
	}

	HipCenter hip(points);
	hip.GetHipCenterBySphere(point_out);
}