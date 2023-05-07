#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImage.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRayCastInterpolateImageFunction.h"

#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkThresholdImageFilter.h"

#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryThresholdImageFilter.h"

#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkConnectedThresholdImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkPointSet.h"
#include "itkImageSeriesWriter.h"
#include "itkImageFileWriter.h"
#include "itkImageDuplicator.h"
#include "itkBinaryFillholeImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkWatershedImageFilter.h"
#include "itkMorphologicalWatershedImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkWhiteTopHatImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkOtsuMultipleThresholdsImageFilter.h"
#include "itkRigid3DTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkNumericTraits.h "
#include "itkNeighborhoodConnectedImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include <opencv2/calib3d/calib3d.hpp>
#include "itkImageToVTKImageFilter.hpp"
#include "vtkFlyingEdges3D.h"

#include "SegmentationException.hpp"
#include "ManualSegmentation.hpp"

using namespace TKA::SEGMENTATION;

const double PI = acos(-1.0);

namespace SegmentationUtils {

    itk::Matrix< double, 3, 3 > getITKTransformFromCVRotation(const cv::Mat& pRotation)
    {
        itk::Matrix< double, 3, 3 > rotation;

        rotation[0][0] = pRotation.at<double>(0, 0);
        rotation[0][1] = pRotation.at<double>(0, 1);
        rotation[0][2] = pRotation.at<double>(0, 2);

        rotation[1][0] = pRotation.at<double>(1, 0);
        rotation[1][1] = pRotation.at<double>(1, 1);
        rotation[1][2] = pRotation.at<double>(1, 2);

        rotation[2][0] = pRotation.at<double>(2, 0);
        rotation[2][1] = pRotation.at<double>(2, 1);
        rotation[2][2] = pRotation.at<double>(2, 2);

        return rotation;
    }

    double getAngleBetweenVectors(const cv::Point3d& a, const cv::Point3d& b)
    {
        double scalar = a.dot(b);
        double magnitude = sqrt((a.dot(a)) * (b.dot(b)));
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

    cv::Mat getRotateMatrix(const cv::Point3d& axis, double angle)
    {
        cv::Mat I = cv::Mat::eye(3, 3, CV_64F);

        if (angle == 0 || angle == PI)
        {
            return I;
        }

        cv::Mat rotationMatrix(3, 3, CV_64F);

        cv::Point3d normaliceAxis = axis;
        normaliceAxis = normaliceAxis / sqrt(normaliceAxis.dot(normaliceAxis));

        rotationMatrix.at <double>(0, 0) = 0;
        rotationMatrix.at <double>(1, 0) = normaliceAxis.z;
        rotationMatrix.at <double>(2, 0) = -normaliceAxis.y;

        rotationMatrix.at <double>(0, 1) = -normaliceAxis.z;
        rotationMatrix.at <double>(1, 1) = 0;
        rotationMatrix.at <double>(2, 1) = normaliceAxis.x;

        rotationMatrix.at <double>(0, 2) = normaliceAxis.y;
        rotationMatrix.at <double>(1, 2) = -normaliceAxis.x;
        rotationMatrix.at <double>(2, 2) = 0;


        cv::Mat result = I + sin(angle)*rotationMatrix + (1.0 - cos(angle))*(rotationMatrix * rotationMatrix);
        return result;
    }

};


ManualSegmentation::ManualSegmentation(const ImageType::Pointer pImage)
{
    mLegImage = pImage;
}

SegmentImageType::Pointer ManualSegmentation::CastImage(const ImageType::Pointer pInput)
{
    constexpr unsigned int Dimension = 3;

    using InputImageType = itk::Image<PixelType, Dimension>;
    using OutputImageType = itk::Image<SegmentPixelType, Dimension>;

    using FilterType = itk::CastImageFilter<InputImageType, OutputImageType>;
    FilterType::Pointer filter = FilterType::New();

    filter->SetInput(pInput);
    filter->Update();
    return filter->GetOutput();
}

ImageType::Pointer ManualSegmentation::GradientAnisotropicDiffusion(const ImageType::Pointer pImage) const
{
    const int numberOfIterations = 5;
    const float timeStep = 0.06;
    const float conductance = 2.0;

    constexpr unsigned int Dimension = 3;
    using OutPixelType = float;
    using OutputImageType = itk::Image<OutPixelType, Dimension>;

    using FilterCast1 = itk::CastImageFilter<ImageType, OutputImageType>;
    FilterCast1::Pointer filterCast1 = FilterCast1::New();

    filterCast1->SetInput(pImage);
    filterCast1->Update();

    using FilterType = itk::GradientAnisotropicDiffusionImageFilter<OutputImageType, OutputImageType>;
    auto filter = FilterType::New();
    filter->SetInput(filterCast1->GetOutput());
    filter->SetNumberOfIterations(numberOfIterations);
    filter->SetTimeStep(timeStep);
    filter->SetConductanceParameter(conductance);
    filter->GlobalWarningDisplayOff();
    filter->Update();

    using FilterCast2 = itk::CastImageFilter<OutputImageType, ImageType>;
    FilterCast2::Pointer filterCast2 = FilterCast2::New();
    filterCast2->SetInput(filter->GetOutput());
    filterCast2->Update();

    return filterCast2->GetOutput();
}

template<typename ImageType>
typename ImageType::Pointer ManualSegmentation::CloneImage(const typename ImageType::Pointer pInput)
{
    using DuplicatorType = itk::ImageDuplicator<ImageType>;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(pInput);
    duplicator->Update();
    return duplicator->GetOutput();
}

template<typename ImageType>
typename ImageType::Pointer ManualSegmentation::FillHole(typename ImageType::Pointer pImage, typename ImageType::PixelType pValue)
{
    using FillholeFilterType = itk::BinaryFillholeImageFilter<ImageType>;
    FillholeFilterType::Pointer fillHoleFilter = FillholeFilterType::New();
    fillHoleFilter->SetInput(pImage);
    fillHoleFilter->SetForegroundValue(pValue);
    fillHoleFilter->Update();
    return fillHoleFilter->GetOutput();
}

void ManualSegmentation::DeepCopy(const ImageType::Pointer pInput, ImageType::Pointer pOutput) const
{
    pOutput->SetRegions(pInput->GetLargestPossibleRegion());
    pOutput->Allocate();

    itk::ImageRegionConstIterator<ImageType> inputIterator(pInput, pInput->GetLargestPossibleRegion());
    itk::ImageRegionIterator<ImageType>      outputIterator(pOutput, pOutput->GetLargestPossibleRegion());

    while (!inputIterator.IsAtEnd())
    {
        outputIterator.Set(inputIterator.Get());
        ++inputIterator;
        ++outputIterator;
    }
}

ImageType::Pointer ManualSegmentation::ApplyThresholds(short pMinthreshold, short pMaxthreshold, bool pDenoiseBefore)
{
    using BinaryThresholdFilterType = itk::BinaryThresholdImageFilter<ImageType, ImageType>;
    BinaryThresholdFilterType::Pointer BinaryFilter = BinaryThresholdFilterType::New();

    if (pDenoiseBefore == true)
    {
        auto tLegImageClean = GradientAnisotropicDiffusion(mLegImage);
        BinaryFilter->SetInput(tLegImageClean);
    }
    else
    {
        BinaryFilter->SetInput(mLegImage);
    }
    
    const ImageType::PixelType  outsidevalue = BinaryBackgroundPixel;
    const ImageType::PixelType  insidevalue = BinaryNoBackgroundPixel;

    BinaryFilter->SetOutsideValue(outsidevalue);
    BinaryFilter->SetInsideValue(insidevalue);

    BinaryFilter->SetUpperThreshold(pMaxthreshold);
    BinaryFilter->SetLowerThreshold(pMinthreshold);
    BinaryFilter->Update();// Running Filter;

    return BinaryFilter->GetOutput();
}


void ManualSegmentation::Get2DSlice(const ImageType::Pointer& pImageIn, ImageType::Pointer& pImageOut, int pSliceNumber, const ViewPlane& pPlane) const
{
    typedef ImageType InputImageType;
    typedef ImageType SliceImageType;
    typedef itk::ExtractImageFilter< InputImageType, SliceImageType > FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetDirectionCollapseToSubmatrix();
    InputImageType::RegionType inputRegion = pImageIn->GetLargestPossibleRegion();
    InputImageType::SizeType size = inputRegion.GetSize();
    size[pPlane] = 1;
    InputImageType::IndexType start = inputRegion.GetIndex();
    start[pPlane] = pSliceNumber;

    InputImageType::RegionType desiredRegion;
    desiredRegion.SetSize(size);
    desiredRegion.SetIndex(start);
    filter->SetExtractionRegion(desiredRegion);
    filter->SetInput(pImageIn);
    filter->Update();
    pImageOut = filter->GetOutput();
}

ImageType::Pointer ManualSegmentation::MakeRotation(const ImageType::Pointer& pImageIn, const ViewPlane& pAnatomicalPlane, const itk::Vector<double, 3>& pObliqueVector, itk::Matrix< double, 3, 3 >& pRotationOut, ImageType::PixelType pDefaultPixelValue)
{
    cv::Point3d currentVector, targetVector;

    if (pAnatomicalPlane == ViewPlane::kSagital_X)
    {
        targetVector = { 1, 0, 0 };
    }
    else if (pAnatomicalPlane == ViewPlane::kCoronal_Y)
    {
        targetVector = { 0, 1, 0 };
    }
    else
    {
        targetVector = { 0, 0, 1 };
    }

    currentVector.x = pObliqueVector[0];
    currentVector.y = pObliqueVector[1];
    currentVector.z = pObliqueVector[2];

    double angle1 = SegmentationUtils::getAngleBetweenVectors(currentVector, targetVector);
    double angle2 = SegmentationUtils::getAngleBetweenVectors(-currentVector, targetVector);
    double angle;

    if (angle1 < angle2)
    {
        angle = angle1;
    }
    else
    {
        angle = angle2;
        currentVector = -currentVector;
    }

    cv::Point3d axis = currentVector.cross(targetVector);

    cv::Mat tempRotation = SegmentationUtils::getRotateMatrix(axis, angle);

    pRotationOut = SegmentationUtils::getITKTransformFromCVRotation(tempRotation);

    /*
    using TransformType = itk::Euler3DTransform<double>;
    TransformType::Pointer initialTransform = TransformType::New();

    ImageType::SpacingType fixedSpacing = pImageIn->GetSpacing();
    itk::Point<ImageType::PixelType, 3>  fixedOrigin = pImageIn->GetOrigin();
    ImageType::RegionType  fixedRegion = pImageIn->GetLargestPossibleRegion();
    ImageType::SizeType fixedSize = fixedRegion.GetSize();

    TransformType::InputPointType centerFixed;
    centerFixed[0] = fixedOrigin[0] + fixedSpacing[0] * fixedSize[0] / 2.0;
    centerFixed[1] = fixedOrigin[1] + fixedSpacing[1] * fixedSize[1] / 2.0;
    centerFixed[2] = fixedOrigin[2] + fixedSpacing[2] * fixedSize[2] / 2.0;

    initialTransform->SetCenter(centerFixed);
    initialTransform->SetMatrix(pRotationOut);

    using ResampleImageFilterType = itk::ResampleImageFilter<ImageType, ImageType>;
    auto resampleFilter = ResampleImageFilterType::New();

    using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, double>;
    auto interpolator = InterpolatorType::New();

    resampleFilter->SetInterpolator(interpolator);
    resampleFilter->SetTransform(initialTransform);
    resampleFilter->SetInput(pImageIn);
    resampleFilter->SetSize(pImageIn->GetLargestPossibleRegion().GetSize());

    resampleFilter->SetOutputOrigin(fixedOrigin);
    resampleFilter->SetOutputSpacing(fixedSpacing);
    resampleFilter->SetDefaultPixelValue(pDefaultPixelValue);

    resampleFilter->Update();
    return resampleFilter->GetOutput();
    */

    ImageType::Pointer result = RotateImage<ImageType>(pImageIn, pRotationOut, pDefaultPixelValue);
    //itk::Matrix< double, 3, 3 > invMatrix = pRotationOut.GetInverse();
    //auto result2 = RotateImage<ImageType>(result, invMatrix, pDefaultPixelValue);
    return result;
}

vtkSmartPointer<vtkPolyData> ManualSegmentation::BinaryImageToPolyData(const SegmentImageType::Pointer pBinaryImage)
{
    using FilterType = itk::ImageToVTKImageFilter<SegmentImageType>;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(pBinaryImage);
    filter->Update();

    auto surface = vtkSmartPointer<vtkFlyingEdges3D>::New();

    surface->SetInputData(filter->GetOutput());

    surface->SetNumberOfContours(1);
    surface->SetValue(0, 1);
    surface->Update();
    return surface->GetOutput();
}

template<typename ImageType>
typename ImageType::Pointer ManualSegmentation::RotateImage(typename ImageType::Pointer pImage, const itk::Matrix< double, 3, 3 >& pRotationIn, typename ImageType::PixelType pDefaultPixelValue)
{
    using TransformType = itk::Euler3DTransform<double>;
    TransformType::Pointer initialTransform = TransformType::New();

    ImageType::SpacingType fixedSpacing = pImage->GetSpacing();
    itk::Point<ImageType::PixelType, 3>  fixedOrigin = pImage->GetOrigin();
    ImageType::RegionType  fixedRegion = pImage->GetLargestPossibleRegion();
    ImageType::SizeType fixedSize = fixedRegion.GetSize();

    TransformType::InputPointType centerFixed;
    centerFixed[0] = fixedOrigin[0] + fixedSpacing[0] * fixedSize[0] / 2.0;
    centerFixed[1] = fixedOrigin[1] + fixedSpacing[1] * fixedSize[1] / 2.0;
    centerFixed[2] = fixedOrigin[2] + fixedSpacing[2] * fixedSize[2] / 2.0;

    initialTransform->SetCenter(centerFixed);
    initialTransform->SetMatrix(pRotationIn);

    using ResampleImageFilterType = itk::ResampleImageFilter<ImageType, ImageType>;
    auto resampleFilter = ResampleImageFilterType::New();

    using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, double>;
    auto interpolator = InterpolatorType::New();

    resampleFilter->SetInterpolator(interpolator);
    resampleFilter->SetTransform(initialTransform);
    resampleFilter->SetInput(pImage);
    resampleFilter->SetSize(pImage->GetLargestPossibleRegion().GetSize());

    resampleFilter->SetOutputOrigin(fixedOrigin);
    resampleFilter->SetOutputSpacing(fixedSpacing);
    resampleFilter->SetDefaultPixelValue(pDefaultPixelValue);

    resampleFilter->Update();
    return resampleFilter->GetOutput();
}

SegmentImageType::Pointer ManualSegmentation::RestoreImageSegmentation(const ImageType::Pointer& pSegmentationIn, const itk::Matrix< double, 3, 3 >& pRotationIn)
{
    ImageType::Pointer fillImage = FillHole<ImageType>(pSegmentationIn, 255);
    itk::Matrix< double, 3, 3 > invMatrix = pRotationIn.GetInverse();
    ImageType::Pointer result = RotateImage<ImageType>(fillImage, invMatrix, 0);
    SegmentImageType::Pointer castImage = CastImage(result);

    return castImage;
}

ImageType::Pointer ManualSegmentation::FloodFillBinaryImageSlice(const ImageType::Pointer& pImage, int8_t pNoBackgroundPixel)
{
    ///////////////////////////// Inverting Background

    using InvertIntensityImageFilterType = itk::InvertIntensityImageFilter<ImageType>;
    auto invertIntensityFilter = InvertIntensityImageFilterType::New();
    invertIntensityFilter->SetInput(pImage);
    invertIntensityFilter->SetMaximum(pNoBackgroundPixel);
    invertIntensityFilter->Update();
    ImageType::Pointer invertImage = invertIntensityFilter->GetOutput();

    //////////////////////////////////////////////////////////////

    const unsigned int Dimension = 3;
    typedef unsigned short                                LabelType;
    typedef itk::Image< LabelType, Dimension >            OutputImageType;
    typedef itk::ShapeLabelObject< LabelType, Dimension > ShapeLabelObjectType;
    typedef itk::LabelMap< ShapeLabelObjectType >         LabelMapType;

    typedef itk::ConnectedComponentImageFilter <ImageType, OutputImageType > ConnectedComponentImageFilterType;
    ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
    connected->SetInput(invertImage);
    connected->Update();
    typedef itk::LabelImageToShapeLabelMapFilter< OutputImageType, LabelMapType> I2LType;
    I2LType::Pointer i2l = I2LType::New();
    i2l->SetInput(connected->GetOutput());
    i2l->Update();
    LabelMapType *labelMap = i2l->GetOutput();

    LabelType background = 0;
    itk::SizeValueType pixels = 0;

    for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n)
    {
        ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);

        if (labelObject->GetNumberOfPixels() > pixels)
        {
            pixels = labelObject->GetNumberOfPixels();
            background = labelObject->GetLabel();
        }
    }

    ImageType::RegionType inputRegion = invertImage->GetLargestPossibleRegion();
    ImageType::SizeType dims = inputRegion.GetSize();
    LabelType* labelImageImBuffer = connected->GetOutput()->GetBufferPointer();

    using DuplicatorType = itk::ImageDuplicator<ImageType>;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(invertImage);
    duplicator->Update();
    ImageType::Pointer result = duplicator->GetOutput();
    ImageType::PixelType* imBuffer = result->GetBufferPointer();

    for (int k = 0; k < dims[2]; k++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int i = 0; i < dims[0]; i++)
            {
                int index = k * dims[1] * dims[0] + j * dims[0] + i;

                if (labelImageImBuffer[index] == background)
                {
                    imBuffer[index] = BinaryBackgroundPixel;
                }
                else
                {
                    imBuffer[index] = BinaryNoBackgroundPixel;
                }
            }
        }
    }
    return result;
}


ImageType2D::Pointer ManualSegmentation::FloodFillBinaryImageSlice(const ImageType2D::Pointer& pImage, int8_t pNoBackgroundPixel)
{
    ///////////////////////////// Inverting Background

    using InvertIntensityImageFilterType = itk::InvertIntensityImageFilter<ImageType2D>;
    auto invertIntensityFilter = InvertIntensityImageFilterType::New();
    invertIntensityFilter->SetInput(pImage);
    invertIntensityFilter->SetMaximum(pNoBackgroundPixel);
    invertIntensityFilter->Update();
    ImageType2D::Pointer invertImage = invertIntensityFilter->GetOutput();

    //////////////////////////////////////////////////////////////

    const unsigned int Dimension = 2;
    typedef unsigned short                                LabelType;
    typedef itk::Image< LabelType, Dimension >            OutputImageType;
    typedef itk::ShapeLabelObject< LabelType, Dimension > ShapeLabelObjectType;
    typedef itk::LabelMap< ShapeLabelObjectType >         LabelMapType;

    typedef itk::ConnectedComponentImageFilter <ImageType2D, OutputImageType > ConnectedComponentImageFilterType;
    ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
    connected->SetInput(invertImage);
    connected->Update();
    typedef itk::LabelImageToShapeLabelMapFilter< OutputImageType, LabelMapType> I2LType;
    I2LType::Pointer i2l = I2LType::New();
    i2l->SetInput(connected->GetOutput());
    i2l->Update();
    LabelMapType *labelMap = i2l->GetOutput();

    LabelType background = 0;
    itk::SizeValueType pixels = 0;

    for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n)
    {
        ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);

        if (labelObject->GetNumberOfPixels() > pixels)
        {
            pixels = labelObject->GetNumberOfPixels();
            background = labelObject->GetLabel();
        }
    }

    ImageType2D::RegionType inputRegion = invertImage->GetLargestPossibleRegion();
    ImageType2D::SizeType dims = inputRegion.GetSize();
    LabelType* labelImageImBuffer = connected->GetOutput()->GetBufferPointer();

    using DuplicatorType = itk::ImageDuplicator<ImageType2D>;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(invertImage);
    duplicator->Update();
    ImageType2D::Pointer result = duplicator->GetOutput();
    ImageType2D::PixelType* imBuffer = result->GetBufferPointer();

    for (int j = 0; j < dims[1]; j++)
    {
        for (int i = 0; i < dims[0]; i++)
        {
            int index = j * dims[0] + i;

            if (labelImageImBuffer[index] == background)
            {
                imBuffer[index] = BinaryBackgroundPixel;
            }
            else
            {
                imBuffer[index] = BinaryNoBackgroundPixel;
            }
        }
    }

    return result;
}
