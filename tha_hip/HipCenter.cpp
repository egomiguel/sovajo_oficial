#include <iostream>
#include <math.h>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <random>
#include "HipCenter.hpp"
#include "HipException.hpp"

using namespace THA::HIP;

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

cv::Mat HipCenter::LeastSquareSolve(const cv::Mat& A, const cv::Mat& B) const
{
	cv::SVD svd(A);
	cv::Mat pinvA = svd.vt.t()* cv::Mat::diag(1. / svd.w)*svd.u.t();
	cv::Mat X = pinvA * cv::Mat(B);
	return X;
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
	cv::Mat X = LeastSquareSolve(A, B);

	point_out[0] = X.at<double>(0, 0) / 2.0;
	point_out[1] = X.at<double>(1, 0) / 2.0;
	point_out[2] = X.at<double>(2, 0) / 2.0;
}

HipCenter::Sphere HipCenter::GetHipCenterBySphere() const
{
	if (is_ellipse_ == true)
	{
        throw HipExceptionCode::SPHERE_HAVE_NOT_BEEN_INITIALIZED;
	}

	HipCenter::Sphere sphere1 = GetSphereCenter(point_list_);
	auto filteredPoints = removeOutliersMAD(point_list_, sphere1.center, sphere1.radius, 3.0);
	HipCenter::Sphere sphere2 = GetSphereCenter(filteredPoints);
	return sphere2;
}

HipCenter::Sphere HipCenter::GetSphereCenter(const std::vector<cv::Point3d>& sphere_points, bool error) const
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
	cv::Mat X = LeastSquareSolve(A, B);

	cv::Point3d center = cv::Point3d(X.at<double>(0, 0) / 2.0, X.at<double>(1, 0) / 2.0, X.at<double>(2, 0) / 2.0);
	double radius =
		std::sqrt(center.x*center.x +
			center.y*center.y +
			center.z*center.z +
			X.at<double>(3, 0));

	HipCenter::Sphere mySphere;
	mySphere.center = center;
	mySphere.radius = radius;

	if (error == false)
	{
		mySphere.error = -1;
	}
	else
	{

	}

	return mySphere;
}

HipCenter::Sphere HipCenter::GetHipCenterFromFile(std::string path)
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
	return hip.GetHipCenterBySphere();
}

std::vector<double> HipCenter::computeResiduals(
	const std::vector<cv::Point3d>& points,
	const cv::Point3d& center,
	double radius)
{
	std::vector<double> residuals;
	residuals.reserve(points.size());

	for (const auto& p : points)
	{
		double dx = p.x - center.x;
		double dy = p.y - center.y;
		double dz = p.z - center.z;

		double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

		residuals.push_back(std::abs(dist - radius));
	}

	return residuals;
}

double HipCenter::computeMAD(const std::vector<double>& residuals)
{
	double med = median(residuals);

	std::vector<double> deviations;
	deviations.reserve(residuals.size());

	for (double r : residuals)
		deviations.push_back(std::abs(r - med));

	return median(deviations);
}

double HipCenter::median(std::vector<double> values)
{
	if (values.empty())
		return 0.0;

	std::sort(values.begin(), values.end());

	size_t n = values.size();

	if (n % 2 == 1)
	{
		return values[n / 2];
	}

	return 0.5 * (values[n / 2 - 1] + values[n / 2]);
}

std::vector<cv::Point3d> HipCenter::removeOutliersMAD(
	const std::vector<cv::Point3d>& points,
	const cv::Point3d& center,
	double radius,
	double zThreshold,
	double maxRemovalFraction)
{
	if (points.empty())
		return {};

	auto residuals = computeResiduals(points, center, radius);

	double med = median(residuals);
	double mad = computeMAD(residuals);

	if (mad < 1e-10)
		return points;

	double sigmaRobust = 1.4826 * mad;

	struct Candidate
	{
		size_t index;
		double absZ;
	};

	std::vector<Candidate> candidates;

	for (size_t i = 0; i < residuals.size(); ++i)
	{
		double z = (residuals[i] - med) / sigmaRobust;

		if (std::abs(z) > zThreshold)
		{
			candidates.push_back(
				{
					i,
					std::abs(z)
				});
		}
	}

	std::sort(
		candidates.begin(),
		candidates.end(),
		[](const Candidate& a, const Candidate& b)
	{
		return a.absZ > b.absZ;
	});

	size_t maxToRemove =
		static_cast<size_t>(
			std::floor(points.size() * maxRemovalFraction));

	if (maxToRemove == 0 && points.size() > 20)
		maxToRemove = 1;

	if (candidates.size() > maxToRemove)
		candidates.resize(maxToRemove);

	std::vector<bool> removeMask(points.size(), false);

	for (const auto& c : candidates)
		removeMask[c.index] = true;

	std::vector<cv::Point3d> filtered;
	filtered.reserve(points.size());

	for (size_t i = 0; i < points.size(); ++i)
	{
		if (!removeMask[i])
			filtered.push_back(points[i]);
	}

	return filtered;
}

HipCenter::Sphere HipCenter::EstimateSphereUncertainty(const std::vector<cv::Point3d>& points, int num_blocks)
{
	std::vector<double> cx, cy, cz, cr;

	size_t N = points.size();
	size_t block_size = N / num_blocks;

	if (block_size < 1)
		block_size = 1;

	for (int b = 0; b < num_blocks; ++b)
	{
		std::vector<cv::Point3d> subset;
		subset.reserve(N - block_size);

		size_t start = b * block_size;
		size_t end = std::min(start + block_size, N);

		for (size_t i = 0; i < N; ++i)
		{
			if (i < start || i >= end)
				subset.push_back(points[i]);
		}

		Sphere s = FitSphere(subset);

		cx.push_back(s.center.x);
		cy.push_back(s.center.y);
		cz.push_back(s.center.z);
		cr.push_back(s.radius);
	}

	// -----------------------------
	// std dev helper
	// -----------------------------
	auto stddev = [](const std::vector<double>& v)
	{
		double mean = 0.0;

		for (double x : v)
			mean += x;

		mean /= v.size();

		double var = 0.0;

		for (double x : v)
			var += (x - mean) * (x - mean);

		var /= (v.size() - 1);

		return std::sqrt(var);
	};

	// -----------------------------
	// resultados
	// -----------------------------
	result.center_std = cv::Point3d(
		stddev(cx),
		stddev(cy),
		stddev(cz));

	result.radius_std = stddev(cr);

	result.center_global_error =
		std::sqrt(
			result.center_std.x * result.center_std.x +
			result.center_std.y * result.center_std.y +
			result.center_std.z * result.center_std.z);

	return result;
}