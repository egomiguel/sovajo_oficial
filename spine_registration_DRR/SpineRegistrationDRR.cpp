#include "SpineRegistrationDRR.hpp"
//#include "RegistrationException.hpp"
//#include "RegistrationPrivate.hpp"
//#include "LeastSquaresICP.hpp"
#include <itkVersorRigid3DTransform.h>
#include <itkImageRegistrationMethodv4.h>
#include <itkMeanSquaresImageToImageMetricv4.h>
#include <itkVersorRigid3DTransform.h>
#include <itkCenteredTransformInitializer.h>
#include <itkRegularStepGradientDescentOptimizerv4.h>
#include <itkResampleImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
//#include <Eigen/Dense>
//#include "itkEuler3DTransform.h"
//#include "itkNormalizedCorrelationTwoImageToOneImageMetric.h"
//#include "itkSiddonJacobsRayCastInterpolateImageFunction.h"
//#include "itkPowellOptimizer.h"
//#include "itkResampleImageFilter.h"
//#include "itkCastImageFilter.h"
//#include "itkRescaleIntensityImageFilter.h"
//#include "itkFlipImageFilter.h"
//#include "itkTwoProjectionImageRegistrationMethod.h"

using namespace SPINE::REGISTRATION_DRR;

SpineRegistrationDRR::SpineRegistrationDRR()
{}

itk::Rigid3DTransform<double>::Pointer SpineRegistrationDRR::ImageRegistration2D3D(const RegistrationImageType::Pointer& pFixedImage, const RegistrationImageType::Pointer& pMovingImage, RegistrationImageType::Pointer& pMovingImageOutput, int pDefaultPixel)
{
	/*
	using TransformType = itk::Euler3DTransform<double>;

	using OptimizerType = itk::PowellOptimizer;

	using MetricType = itk::NormalizedCorrelationTwoImageToOneImageMetric<RegistrationDRRInternalPixelType, RegistrationDRRInternalPixelType>;

	using InterpolatorType = itk::SiddonJacobsRayCastInterpolateImageFunction<RegistrationDRRInternalPixelType, double>;

	using RegistrationType = itk::TwoProjectionImageRegistrationMethod<RegistrationDRRInternalPixelType, RegistrationDRRInternalPixelType>;

	MetricType::Pointer       metric = MetricType::New();
	TransformType::Pointer    transform = TransformType::New();
	OptimizerType::Pointer    optimizer = OptimizerType::New();
	InterpolatorType::Pointer interpolator1 = InterpolatorType::New();
	InterpolatorType::Pointer interpolator2 = InterpolatorType::New();
	RegistrationType::Pointer registration = RegistrationType::New();

	metric->ComputeGradientOff();
	metric->SetSubtractMean(true);

	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);
	registration->SetTransform(transform);
	registration->SetInterpolator1(interpolator1);
	registration->SetInterpolator2(interpolator2);


	// To simply Siddon-Jacob's fast ray-tracing algorithm, we force the origin of the CT image
	// to be (0,0,0). Because we align the CT isocenter with the central axis, the projection
	// geometry is fully defined. The origin of the CT image becomes irrelavent.
	RegistrationDRRImageType3D::Pointer image3DIn = pImage3D;
	RegistrationDRRImageType3D::PointType image3DOrigin;
	image3DOrigin[0] = 0.0;
	image3DOrigin[1] = 0.0;
	image3DOrigin[2] = 0.0;
	image3DIn->SetOrigin(image3DOrigin);

	// The input 2D images were loaded as 3D images. They were considered
	// as a single slice from a 3D volume. By default, images stored on the
	// disk are treated as if they have RAI orientation. After view point
	// transformation, the order of 2D image pixel reading is equivalent to
	// from inferior to superior. This is contradictory to the traditional
	// 2D x-ray image storage, in which a typical 2D image reader reads and
	// writes images from superior to inferior. Thus the loaded 2D DICOM
	// images should be flipped in y-direction. This was done by using a.
	// FilpImageFilter.
	using FlipFilterType = itk::FlipImageFilter<RegistrationDRRInternalPixelType>;
	FlipFilterType::Pointer flipFilter1 = FlipFilterType::New();
	//FlipFilterType::Pointer flipFilter2 = FlipFilterType::New();

	using FlipAxesArrayType = FlipFilterType::FlipAxesArrayType;
	FlipAxesArrayType flipArray;
	flipArray[0] = false;
	flipArray[1] = true;
	flipArray[2] = false;

	flipFilter1->SetFlipAxes(flipArray);
	//flipFilter2->SetFlipAxes(flipArray);

	flipFilter1->SetInput(pImage2D);
	//flipFilter2->SetInput(imageReader2D2->GetOutput());

	// The input 2D images may have 16 bits. We rescale the pixel value to between 0-255.
	using Input2DRescaleFilterType = itk::RescaleIntensityImageFilter<RegistrationDRRInternalPixelType, RegistrationDRRInternalPixelType>;

	Input2DRescaleFilterType::Pointer rescaler2D1 = Input2DRescaleFilterType::New();
	rescaler2D1->SetOutputMinimum(0);
	rescaler2D1->SetOutputMaximum(255);
	rescaler2D1->SetInput(flipFilter1->GetOutput());


	//  The 3D CT dataset is casted to the internal image type using
	//  {CastImageFilters}.

	using CastFilterType3D = itk::CastImageFilter<RegistrationDRRImageType3D, RegistrationDRRInternalPixelType>;

	CastFilterType3D::Pointer caster3D = CastFilterType3D::New();
	caster3D->SetInput(image3DIn);

	rescaler2D1->Update();
	//rescaler2D2->Update();
	caster3D->Update();


	registration->SetFixedImage1(rescaler2D1->GetOutput());
	//registration->SetFixedImage2(rescaler2D2->GetOutput());
	registration->SetMovingImage(caster3D->GetOutput());

	// Initialise the transform
	// ~~~~~~~~~~~~~~~~~~~~~~~~

	// Set the order of the computation. Default ZXY
	transform->SetComputeZYX(true);


	// The transform is initialised with the translation [tx,ty,tz] and
	// rotation [rx,ry,rz] specified on the command line

	TransformType::OutputVectorType translation;

	translation[0] = 0;
	translation[1] = 0;
	translation[2] = 0;

	transform->SetTranslation(translation);

	// constant for converting degrees to radians
	const double dtr = (atan(1.0) * 4.0) / 180.0;
	transform->SetRotation(dtr * 0, dtr * 0, dtr * 0);

	// The centre of rotation is set by default to the centre of the 3D
	// volume but can be offset from this position using a command
	// line specified translation [cx,cy,cz]

	RegistrationDRRImageType3D::PointType origin3D = image3DIn->GetOrigin();
	const itk::Vector<double, 3> resolution3D = image3DIn->GetSpacing();

	using ImageRegionType3D = RegistrationDRRImageType3D::RegionType;
	using SizeType3D = ImageRegionType3D::SizeType;

	ImageRegionType3D region3D = caster3D->GetOutput()->GetBufferedRegion();
	SizeType3D        size3D = region3D.GetSize();

	TransformType::InputPointType isocenter;
	// Set the center of the image as the isocenter.
	isocenter[0] = origin3D[0] + resolution3D[0] * static_cast<double>(size3D[0]) / 2.0;
	isocenter[1] = origin3D[1] + resolution3D[1] * static_cast<double>(size3D[1]) / 2.0;
	isocenter[2] = origin3D[2] + resolution3D[2] * static_cast<double>(size3D[2]) / 2.0;
	
	transform->SetCenter(isocenter);

	if (true)
	{
		std::cout << "3D image size: " << size3D[0] << ", " << size3D[1] << ", " << size3D[2] << std::endl
			<< "   resolution: " << resolution3D[0] << ", " << resolution3D[1] << ", " << resolution3D[2] << std::endl
			<< "Transform: " << transform << std::endl;
	}


	// Set the origin of the 2D image
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// For correct (perspective) projection of the 3D volume, the 2D
	// image needs to be placed at a certain distance (the source-to-
	// isocenter distance {scd} ) from the focal point, and the normal
	// from the imaging plane to the focal point needs to be specified.
	//
	// By default, the imaging plane normal is set by default to the
	// center of the 2D image but may be modified from this using the
	// command line parameters [image1centerX, image1centerY,
	// image2centerX, image2centerY].

	double origin2D1[3];

	// Note: Two 2D images may have different image sizes and pixel dimensions, although
	// scd are the same.

	const itk::Vector<double, 3> resolution2D1 = pImage2D->GetSpacing();

	using ImageRegionType2D = RegistrationDRRInternalPixelType::RegionType;
	using SizeType2D = ImageRegionType2D::SizeType;

	ImageRegionType2D region2D1 = rescaler2D1->GetOutput()->GetBufferedRegion();
	SizeType2D        size2D1 = region2D1.GetSize();

	double image1centerX = 0.0;
	double image1centerY = 0.0;
	double scd = 1000.;

	image1centerX = ((double)size2D1[0] - 1.) / 2.;
	image1centerY = ((double)size2D1[1] - 1.) / 2.;
	

	// 2D Image 1
	origin2D1[0] = -resolution2D1[0] * image1centerX;
	origin2D1[1] = -resolution2D1[1] * image1centerY;
	origin2D1[2] = -scd;

	rescaler2D1->GetOutput()->SetOrigin(origin2D1);

	registration->SetFixedImageRegion1(rescaler2D1->GetOutput()->GetBufferedRegion());


	// Initialize the ray cast interpolator
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// The ray cast interpolator is used to project the 3D volume. It
	// does this by casting rays from the (transformed) focal point to
	// each (transformed) pixel coordinate in the 2D image.
	//
	// In addition a threshold may be specified to ensure that only
	// intensities greater than a given value contribute to the
	// projected volume. This can be used, for instance, to remove soft
	// tissue from projections of CT data and force the registration
	// to find a match which aligns bony structures in the images.

	// 2D Image 1
	double threshold = 0.0;
	interpolator1->SetProjectionAngle(dtr * 1);
	interpolator1->SetFocalPointToIsocenterDistance(scd);
	interpolator1->SetThreshold(threshold);
	interpolator1->SetTransform(transform);

	interpolator1->Initialize();


	// Set up the transform and start position
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// The registration start position is intialised using the
	// transformation parameters.

	registration->SetInitialTransformParameters(transform->GetParameters());

	// We wish to minimize the negative normalized correlation similarity measure.

	// optimizer->SetMaximize( true );  // for GradientDifferenceTwoImageToOneImageMetric
	optimizer->SetMaximize(false); // for NCC

	optimizer->SetMaximumIteration(10);
	optimizer->SetMaximumLineIteration(4); // for Powell's method
	optimizer->SetStepLength(4.0);
	optimizer->SetStepTolerance(0.02);
	optimizer->SetValueTolerance(0.001);

	// The optimizer weightings are set such that one degree equates to
	// one millimeter.

	itk::Optimizer::ScalesType weightings(transform->GetNumberOfParameters());

	weightings[0] = 1. / dtr;
	weightings[1] = 1. / dtr;
	weightings[2] = 1. / dtr;
	weightings[3] = 1.;
	weightings[4] = 1.;
	weightings[5] = 1.;

	optimizer->SetScales(weightings);

	using ParametersType = RegistrationType::ParametersType;
	ParametersType finalParameters = registration->GetLastTransformParameters();

	const double RotationAlongX = finalParameters[0] / dtr; // Convert radian to degree
	const double RotationAlongY = finalParameters[1] / dtr;
	const double RotationAlongZ = finalParameters[2] / dtr;
	const double TranslationAlongX = finalParameters[3];
	const double TranslationAlongY = finalParameters[4];
	const double TranslationAlongZ = finalParameters[5];

	const int numberOfIterations = optimizer->GetCurrentIteration();

	const double bestValue = optimizer->GetValue();

	std::cout << "Result = " << std::endl;
	std::cout << " Rotation Along X = " << RotationAlongX << " deg" << std::endl;
	std::cout << " Rotation Along Y = " << RotationAlongY << " deg" << std::endl;
	std::cout << " Rotation Along Z = " << RotationAlongZ << " deg" << std::endl;
	std::cout << " Translation X = " << TranslationAlongX << " mm" << std::endl;
	std::cout << " Translation Y = " << TranslationAlongY << " mm" << std::endl;
	std::cout << " Translation Z = " << TranslationAlongZ << " mm" << std::endl;
	std::cout << " Number Of Iterations = " << numberOfIterations << std::endl;
	std::cout << " Metric value  = " << bestValue << std::endl;

	TransformType::Pointer finalTransform = TransformType::New();
	// The transform is determined by the parameters and the rotation center.
	finalTransform->SetParameters(finalParameters);
	finalTransform->SetCenter(isocenter);

	using ResampleFilterType = itk::ResampleImageFilter<RegistrationDRRInternalPixelType, RegistrationDRRInternalPixelType>;

	// The ResampleImageFilter is the driving force for the projection image generation.
	ResampleFilterType::Pointer resampleFilter1 = ResampleFilterType::New();

	resampleFilter1->SetInput(caster3D->GetOutput()); // Link the 3D volume.
	resampleFilter1->SetDefaultPixelValue(0);

	// The parameters of interpolator1, such as ProjectionAngle and FocalPointToIsocenterDistance
	// have been set before registration. Here we only need to replace the initial
	// transform with the final transform.
	interpolator1->SetTransform(finalTransform);
	interpolator1->Initialize();
	resampleFilter1->SetInterpolator(interpolator1);

	// The output 2D projection image has the same image size, origin, and the pixel spacing as
	// those of the input 2D image.
	resampleFilter1->SetSize(rescaler2D1->GetOutput()->GetLargestPossibleRegion().GetSize());
	resampleFilter1->SetOutputOrigin(rescaler2D1->GetOutput()->GetOrigin());
	resampleFilter1->SetOutputSpacing(rescaler2D1->GetOutput()->GetSpacing());

	// As explained before, the computed projection is upsided-down.
	// Here we use a FilpImageFilter to flip the images in y-direction.
	flipFilter1->SetInput(resampleFilter1->GetOutput());

	// Rescale the intensity of the projection images to 0-255 for output.
	using RescaleFilterType = itk::RescaleIntensityImageFilter<RegistrationDRRInternalPixelType, RegistrationDRRImageType3D>;

	RescaleFilterType::Pointer rescaler1 = RescaleFilterType::New();
	rescaler1->SetOutputMinimum(0);
	rescaler1->SetOutputMaximum(255);
	rescaler1->SetInput(flipFilter1->GetOutput());

	pImageOutput = rescaler1->GetOutput();

	TransformType::MatrixType matrix = finalTransform->GetMatrix();
	TransformType::OffsetType offset = finalTransform->GetOffset();

	itk::Rigid3DTransform<double>::Pointer result = itk::VersorRigid3DTransform<double>::New();
	result->SetMatrix(matrix);
	result->SetOffset(offset);
	return result;
	*/

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
typename ImageTypeOutput::Pointer SpineRegistrationDRR::castImage(const typename ImageTypeInput::Pointer input)
{
    using FilterType = itk::CastImageFilter<ImageTypeInput, ImageTypeOutput>;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(input);
    filter->Update();
    return filter->GetOutput();
}
