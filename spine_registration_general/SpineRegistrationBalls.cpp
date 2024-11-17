#include "SpineRegistrationBalls.hpp"
#include "RegistrationException.hpp"
#include "LeastSquaresICP.hpp"
#include <itkVersorRigid3DTransform.h>
#include <itkImageRegistrationMethodv4.h>
#include <itkMeanSquaresImageToImageMetricv4.h>
#include <itkVersorRigid3DTransform.h>
#include <itkCenteredTransformInitializer.h>
#include <itkRegularStepGradientDescentOptimizerv4.h>
#include <itkResampleImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <Eigen/Dense>

using namespace SPINE::REGISTRATION;

SpineRegistrationBalls::SpineRegistrationBalls(const std::vector<PointTypeITK>& pTargetBallsOnCT, const std::vector<PointTypeITK>& pSourceExternalBalls)
{
	mTransform = cv::Mat(6, 1, CV_64F);

    if (pTargetBallsOnCT.size() < 4 || pSourceExternalBalls.size() < 4)
    {
        mTransform.at<double>(0, 0) = 0;
        mTransform.at<double>(1, 0) = 0;
        mTransform.at<double>(2, 0) = 0;

        mTransform.at<double>(3, 0) = 0;
        mTransform.at<double>(4, 0) = 0;
        mTransform.at<double>(5, 0) = 0;

        return; //Exception ****************************************************
    }

	cv::Mat mTargetOnCT(pTargetBallsOnCT.size(), 3, CV_64F);
	mSourceExternalBalls = pSourceExternalBalls;

	for (int i = 0; i < pTargetBallsOnCT.size(); i++)
	{
		cv::Point3d newPoint = cv::Point3d(pTargetBallsOnCT[i][0], pTargetBallsOnCT[i][1], pTargetBallsOnCT[i][2]);

		mTargetOnCT.at<double>(i, 0) = newPoint.x;
		mTargetOnCT.at<double>(i, 1) = newPoint.y;
		mTargetOnCT.at<double>(i, 2) = newPoint.z;
	}

	SpineRegistrationBalls::PointsFit targetFit = getBestFit(pTargetBallsOnCT);
	SpineRegistrationBalls::PointsFit sourceFit = getBestFit(pSourceExternalBalls);

	cv::Mat transform = GetGeneralRotateTransformVectors(sourceFit.normal, targetFit.normal);

    ///////////////////////////////////////////////////////////////////

    Eigen::Matrix3d rotation(3, 3);

    rotation(0, 0) = transform.at<double>(0, 0);
    rotation(1, 0) = transform.at<double>(1, 0);
    rotation(2, 0) = transform.at<double>(2, 0);

    rotation(0, 1) = transform.at<double>(0, 1);
    rotation(1, 1) = transform.at<double>(1, 1);
    rotation(2, 1) = transform.at<double>(2, 1);

    rotation(0, 2) = transform.at<double>(0, 2);
    rotation(1, 2) = transform.at<double>(1, 2);
    rotation(2, 2) = transform.at<double>(2, 2);

    Eigen::Vector3d rot = rotation.eulerAngles(0, 1, 2);

	cv::Point3d translation = targetFit.center - sourceFit.center;

    mTransform.at<double>(0, 0) = translation.x;
    mTransform.at<double>(1, 0) = translation.y;
    mTransform.at<double>(2, 0) = translation.z;

    mTransform.at<double>(3, 0) = rot(0);
    mTransform.at<double>(4, 0) = rot(1);
    mTransform.at<double>(5, 0) = rot(2);
}

SpineRegistrationBalls::PointsFit SpineRegistrationBalls::getBestFit(const std::vector<PointTypeITK>& pPoints) const
{
	int tSize = pPoints.size();
	cv::Point3d averagePoint(0, 0, 0);
	std::vector<cv::Point3d> myPoints;

	SpineRegistrationBalls::PointsFit myFit;
	myFit.center = averagePoint;
	myFit.normal = averagePoint;

	for (int i = 0; i < tSize; i++)
	{
		cv::Point3d temp = cv::Point3d(pPoints[i][0], pPoints[i][1], pPoints[i][2]);
		averagePoint = averagePoint + temp;
		myPoints.push_back(temp);
	}

	if (tSize > 0)
	{
		averagePoint = averagePoint / double(tSize);
	}
	else
	{
		return myFit;
	}

	std::vector<cv::Point3d> center(tSize, averagePoint);

	cv::Mat A = cv::Mat(tSize, 3, CV_64F, myPoints.data());
	cv::Mat C(tSize, 3, CV_64F, center.data());

	//Subtracting centroid_ from the set of points to center on the origin.
	cv::Mat CA = A - C;
	cv::SVD svd(CA);
	cv::Mat Last_V_Row = svd.vt.row(2);
	double a = Last_V_Row.at<double>(0);
	double b = Last_V_Row.at<double>(1);
	double c = Last_V_Row.at<double>(2);

	myFit.center = averagePoint;
	myFit.normal = cv::Point3d(a, b, c);

	return myFit;
}

itk::Rigid3DTransform<double>::Pointer SpineRegistrationBalls::GetAlignment(double& pError)
{
    LeastSquaresICP myICP(mSourceExternalBalls);

    pError = myICP.LeastSquaresRandomInit(mTargetOnCT, mTransform);

    Eigen::AngleAxisd init_rotationZ(mTransform.at<double>(5, 0), Eigen::Vector3d::UnitZ());
    Eigen::AngleAxisd init_rotationY(mTransform.at<double>(4, 0), Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd init_rotationX(mTransform.at<double>(3, 0), Eigen::Vector3d::UnitX());

    Eigen::Translation3d init_translation(mTransform.at<double>(0, 0), mTransform.at<double>(1, 0), mTransform.at<double>(2, 0));

    Eigen::Matrix4d matrix = (init_translation * init_rotationX * init_rotationY * init_rotationZ).matrix();

    itk::Matrix< double, 3, 3 > rotation;
    itk::Vector< double, 3 > translate;

    rotation[0][0] = matrix(0, 0);
    rotation[0][1] = matrix(0, 1);
    rotation[0][2] = matrix(0, 2);

    rotation[1][0] = matrix(1, 0);
    rotation[1][1] = matrix(1, 1);
    rotation[1][2] = matrix(1, 2);

    rotation[2][0] = matrix(2, 0);
    rotation[2][1] = matrix(2, 1);
    rotation[2][2] = matrix(2, 2);

    translate[0] = matrix(0, 3);
    translate[1] = matrix(1, 3);
    translate[2] = matrix(2, 3);

    itk::Rigid3DTransform<double>::Pointer result = itk::VersorRigid3DTransform<double>::New();
    result->SetMatrix(rotation);
    result->SetOffset(translate);

    return result;
}

cv::Mat SpineRegistrationBalls::getRotateMatrix(const cv::Point3d& axis, double angle)
{
	cv::Mat rotationMatrix(3, 3, CV_64F);

	cv::Point3d normaliceAxis = axis;
	double norm = sqrt(normaliceAxis.dot(normaliceAxis));
	normaliceAxis = normaliceAxis / norm;

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

cv::Mat SpineRegistrationBalls::GetGeneralRotateTransformVectors(const cv::Point3d pFromVector, const cv::Point3d pToVector)
{
	cv::Point3d fromVector = pFromVector;
	cv::Point3d toVector = pToVector;
	fromVector = fromVector / sqrt(fromVector.dot(fromVector));
	toVector = toVector / sqrt(toVector.dot(toVector));

	cv::Point3d rotationAxis = fromVector.cross(toVector);
	rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));

	double rotationAngle = getAngleBetweenVectors(fromVector, toVector);

	cv::Mat rotate = getRotateMatrix(rotationAxis, rotationAngle);

	return rotate;
}

double SpineRegistrationBalls::getAngleBetweenVectors(const cv::Point3d& a, const cv::Point3d& b)
{
	double scalar = a.dot(b);
	double magnitude = sqrt((a.dot(a)) * (b.dot(b)));
	double tCos = scalar / magnitude;
	if (tCos <= -1.0)
	{
		return acos(-1.0);
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