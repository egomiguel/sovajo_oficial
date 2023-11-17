#include "GeneralRegistrationPointsToPolydata.hpp"
#include "RegistrationException.hpp"
#include "RegistrationPrivate.hpp"
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

using namespace PKA::REGISTRATION;

GeneralRegistrationPointsToPolydata::GeneralRegistrationPointsToPolydata(const vtkSmartPointer<vtkPolyData>& pPoly, const std::vector<PointTypeITK>& pTargetPointsOnPoly, const std::vector<PointTypeITK>& pSourceExternalPoints)
{
    mPoly = pPoly;
	mTransform = cv::Mat(6, 1, CV_64F);

    if (pSourceExternalPoints.size() != pTargetPointsOnPoly.size() || pSourceExternalPoints.size() == 0)
    {
        mTransform.at<double>(0, 0) = 0;
        mTransform.at<double>(1, 0) = 0;
        mTransform.at<double>(2, 0) = 0;

        mTransform.at<double>(3, 0) = 0;
        mTransform.at<double>(4, 0) = 0;
        mTransform.at<double>(5, 0) = 0;

        return;
    }

    Eigen::Matrix4d transform;
    int tSize = pSourceExternalPoints.size();

    ////////////////////////////////////////////////////////////// Find Centroid
    Eigen::Matrix<double, 4, 1> centroid_src, centroid_tgt;

    centroid_src.setZero();
    centroid_tgt.setZero();

    auto it1 = pSourceExternalPoints.begin();
    auto it2 = pSourceExternalPoints.end();

    for (; it1 != it2; ++it1)
    {
        centroid_src[0] += (*it1)[0];
        centroid_src[1] += (*it1)[1];
        centroid_src[2] += (*it1)[2];
    }

    centroid_src /= static_cast<double>(tSize);
    centroid_src[3] = 1;

    it1 = pTargetPointsOnPoly.begin();
    it2 = pTargetPointsOnPoly.end();

    for (; it1 != it2; ++it1)
    {
        centroid_tgt[0] += (*it1)[0];
        centroid_tgt[1] += (*it1)[1];
        centroid_tgt[2] += (*it1)[2];
    }

    centroid_tgt /= static_cast<double>(tSize);
    centroid_tgt[3] = 1;

    //////////////////////////////////////////////////////////////////////////

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> cloud_src_demean, cloud_tgt_demean;
    cloud_src_demean = Eigen::Matrix<double, 4, Eigen::Dynamic>::Zero(4, tSize);
    cloud_tgt_demean = Eigen::Matrix<double, 4, Eigen::Dynamic>::Zero(4, tSize);

    it1 = pSourceExternalPoints.begin();
    it2 = pSourceExternalPoints.end();
    int i = 0;
    while (it1 != it2)
    {
        cloud_src_demean(0, i) = (*it1)[0] - centroid_src[0];
        cloud_src_demean(1, i) = (*it1)[1] - centroid_src[1];
        cloud_src_demean(2, i) = (*it1)[2] - centroid_src[2];
        ++i;
        ++it1;
    }

    it1 = pTargetPointsOnPoly.begin();
    it2 = pTargetPointsOnPoly.end();
    i = 0;
    while (it1 != it2)
    {
        cloud_tgt_demean(0, i) = (*it1)[0] - centroid_tgt[0];
        cloud_tgt_demean(1, i) = (*it1)[1] - centroid_tgt[1];
        cloud_tgt_demean(2, i) = (*it1)[2] - centroid_tgt[2];
        ++i;
        ++it1;
    }

    ////////////////////////////////////////////////////////////////////

    transform = Eigen::Matrix4d::Identity();

    // Assemble the correlation matrix H = source * target'
    Eigen::Matrix<double, 3, 3> H = (cloud_src_demean * cloud_tgt_demean.transpose()).topLeftCorner(3, 3);

    // Compute the Singular Value Decomposition
    Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3> > svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix<double, 3, 3> u = svd.matrixU();
    Eigen::Matrix<double, 3, 3> v = svd.matrixV();

    // Compute R = V * U'
    if (u.determinant() * v.determinant() < 0)
    {
        for (int x = 0; x < 3; ++x)
            v(x, 2) *= -1;
    }

    Eigen::Matrix<double, 3, 3> R = v * u.transpose();

    // Return the correct transformation
    transform.topLeftCorner(3, 3) = R;
    const Eigen::Matrix<double, 3, 1> Rc(R * centroid_src.head(3));
    transform.block(0, 3, 3, 1) = centroid_tgt.head(3) - Rc;

    ///////////////////////////////////////////////////////////////////

    Eigen::Matrix3d rotation(3, 3);

    rotation(0, 0) = transform(0, 0);
    rotation(1, 0) = transform(1, 0);
    rotation(2, 0) = transform(2, 0);

    rotation(0, 1) = transform(0, 1);
    rotation(1, 1) = transform(1, 1);
    rotation(2, 1) = transform(2, 1);

    rotation(0, 2) = transform(0, 2);
    rotation(1, 2) = transform(1, 2);
    rotation(2, 2) = transform(2, 2);

    Eigen::Vector3d rot = rotation.eulerAngles(0, 1, 2);

    mTransform.at<double>(0, 0) = transform(0, 3);
    mTransform.at<double>(1, 0) = transform(1, 3);
    mTransform.at<double>(2, 0) = transform(2, 3);

    mTransform.at<double>(3, 0) = rot(0);
    mTransform.at<double>(4, 0) = rot(1);
    mTransform.at<double>(5, 0) = rot(2);

}

itk::Rigid3DTransform<double>::Pointer GeneralRegistrationPointsToPolydata::MakeFinalAlignment(const std::vector<PointTypeITK>& pAlignmentPoints, double& pError)
{
    LeastSquaresICP myICP(pAlignmentPoints);

    pError = myICP.LeastSquares(mPoly, mTransform);

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

itk::Rigid3DTransform<double>::Pointer GeneralRegistrationPointsToPolydata::MultiResImageRegistration(const RegistrationImageType::Pointer& pFixedImage, const RegistrationImageType::Pointer& pMovingImage, RegistrationImageType::Pointer& pMovingImageOutput, int pDefaultPixel)
{
    using FloatImageType = itk::Image<float, 3>;

    FloatImageType::Pointer fixedImage = castImage<RegistrationImageType, FloatImageType>(pFixedImage);
    FloatImageType::Pointer movingImage = castImage<RegistrationImageType, FloatImageType>(pMovingImage);

    using TransformType = itk::VersorRigid3DTransform<double>;

    using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
    using MetricType = itk::MeanSquaresImageToImageMetricv4<FloatImageType, FloatImageType>;
    using RegistrationType = itk::ImageRegistrationMethodv4<FloatImageType, FloatImageType, TransformType>;

    auto metric = MetricType::New();
    auto optimizer = OptimizerType::New();
    auto registration = RegistrationType::New();

    registration->SetMetric(metric);
    registration->SetOptimizer(optimizer);

    auto initialTransform = TransformType::New();

    registration->SetFixedImage(fixedImage);
    registration->SetMovingImage(movingImage);

    using TransformInitializerType = itk::CenteredTransformInitializer<TransformType, FloatImageType, FloatImageType>;
    auto initializer = TransformInitializerType::New();

    initializer->SetTransform(initialTransform);
    initializer->SetFixedImage(fixedImage);
    initializer->SetMovingImage(movingImage);

    initializer->MomentsOn();

    initializer->InitializeTransform();

    using VersorType = TransformType::VersorType;
    using VectorType = VersorType::VectorType;
    VersorType rotation;
    VectorType axis;
    axis[0] = 0.0;
    axis[1] = 0.0;
    axis[2] = 1.0;
    double angle = 0;
    rotation.Set(axis, angle);
    initialTransform->SetRotation(rotation);

    registration->SetInitialTransform(initialTransform);

    using OptimizerScalesType = OptimizerType::ScalesType;
    OptimizerScalesType optimizerScales(initialTransform->GetNumberOfParameters());
    const double anglesScale = 1.0 / 1000.0; // need to be small
    optimizerScales[0] = 1.0;
    optimizerScales[1] = 1.0;
    optimizerScales[2] = 1.0;
    optimizerScales[3] = anglesScale;
    optimizerScales[4] = anglesScale;
    optimizerScales[5] = anglesScale;
    optimizer->SetScales(optimizerScales);
    optimizer->SetNumberOfIterations(200);
    optimizer->SetLearningRate(0.5);
    optimizer->SetMinimumStepLength(0.01);
    optimizer->SetReturnBestParametersAndValue(true);

    int numberOfLevels = 3;

    RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
    shrinkFactorsPerLevel.SetSize(3);
    shrinkFactorsPerLevel[0] = 3;
    shrinkFactorsPerLevel[1] = 2;
    shrinkFactorsPerLevel[2] = 1;

    RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
    smoothingSigmasPerLevel.SetSize(3);
    smoothingSigmasPerLevel[0] = 0.2;
    smoothingSigmasPerLevel[1] = 0.2;
    smoothingSigmasPerLevel[2] = 0.2;

    registration->SetNumberOfLevels(numberOfLevels);
    registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
    registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);
    registration->Update();
    const TransformType::ParametersType finalParameters = registration->GetOutput()->Get()->GetParameters();

    auto finalTransform = TransformType::New();

    finalTransform->SetFixedParameters(registration->GetOutput()->Get()->GetFixedParameters());
    finalTransform->SetParameters(finalParameters);

    TransformType::MatrixType matrix = finalTransform->GetMatrix();
    TransformType::OffsetType offset = finalTransform->GetOffset();

    using ResampleFilterType = itk::ResampleImageFilter<FloatImageType, FloatImageType>;

    auto resampler = ResampleFilterType::New();

    resampler->SetTransform(finalTransform);
    resampler->SetInput(movingImage);

    resampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
    resampler->SetOutputOrigin(fixedImage->GetOrigin());
    resampler->SetOutputSpacing(fixedImage->GetSpacing());
    resampler->SetOutputDirection(fixedImage->GetDirection());
    resampler->SetDefaultPixelValue(pDefaultPixel);
    resampler->Update();
    auto finalImage = resampler->GetOutput();

    pMovingImageOutput = castImage<FloatImageType, RegistrationImageType>(finalImage);

    itk::Rigid3DTransform<double>::Pointer result = itk::VersorRigid3DTransform<double>::New();
    result->SetMatrix(matrix);
    result->SetOffset(offset);
    return result;
}

template<typename ImageTypeInput, typename ImageTypeOutput>
typename ImageTypeOutput::Pointer GeneralRegistrationPointsToPolydata::castImage(const typename ImageTypeInput::Pointer input)
{
    using FilterType = itk::CastImageFilter<ImageTypeInput, ImageTypeOutput>;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(input);
    filter->Update();
    return filter->GetOutput();
}
