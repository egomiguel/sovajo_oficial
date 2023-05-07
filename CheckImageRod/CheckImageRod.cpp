#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImage.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkCenteredEuler3DTransform.h"
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

#include "itkOtsuMultipleThresholdsImageFilter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkCastImageFilter.h"
#include "itkExtractImageFilter.h"

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

#include "CheckImageRodException.hpp"
#include "CheckImageRod.hpp"

using namespace IMAGE_ROD;

template<typename ImageType>
void SaveImageTest(typename ImageType::Pointer img, std::string name, bool verbose = true)
{
    std::string fullName = name + ".nrrd";
    using WriterType = itk::ImageFileWriter<ImageType>;
    typename WriterType::Pointer writer = WriterType::New();

    writer->SetFileName(fullName);
    writer->SetInput(img);

    if (verbose == true)
    {
        std::cout << "Writing the image" <<std::endl;
    }

    try
    {
        writer->Update();
        if (verbose == true)
        {
            std::cout << "Write image!!!" <<std::endl;
        }
    }
    catch (const itk::ExceptionObject & ex)
    {
        std::cout << ex.what() << std::endl;
    }
}

template<typename ImageType>
typename ImageType::Pointer CloneImageTest(const typename ImageType::Pointer input)
{
    using DuplicatorType = itk::ImageDuplicator<ImageType>;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(input);
    duplicator->Update();
    return duplicator->GetOutput();
}

struct LabelsInfo
{
    unsigned short MaxVolLabel;
    double elongation;
    ImageType::RegionType region;
    LabelsInfo(unsigned short MaxVolLabel, double elong, ImageType::RegionType region)
    {
        this->MaxVolLabel = MaxVolLabel;
        this->region = region;
        ImageType::SizeType tSize = this->region.GetSize();

        double h, base1, base2;

        base1 = tSize[0];
        base2 = tSize[1];
        h = tSize[2];

        if (base1 < 1)
        {
            base1 = 1;
        }

        if (base2 < 1)
        {
            base2 = 1;
        }

        double coef;

        if (base1 > base2)
        {
            coef = h / base1;

            if (coef < 2 && coef >= 1)
            {
                coef = 4;
            }
            else if (coef < 1)
            {
                coef = 8;
            }
            else
            {
                coef = 1;
            }

            this->elongation = h / (coef * (base1 / base2));
        }
        else
        {
            coef = h / base2;

            if (coef < 2 && coef >= 1)
            {
                coef = 4;
            }
            else if (coef < 1)
            {
                coef = 8;
            }
            else
            {
                coef = 1;
            }

            this->elongation = h / (coef * (base2 / base1));
        }
    }

    bool Similarity(ImageType::SizeType pSize)
    {
        ImageType::SizeType tSize = this->region.GetSize();

        int a = (pSize[0] - tSize[0]);
        int b = (pSize[1] - tSize[1]);
        int c = (pSize[2] - tSize[2]);

        if (abs(a) == 0 || abs(b) == 0 || abs(c) == 0)
        {
            return true;
        }
        return false;
    }
};

struct SortByElongation
{
    bool operator()(const LabelsInfo & label1, const LabelsInfo & label2) { return label1.elongation > label2.elongation; }
};

CheckImageRod::CheckImageRod(const ImageType::Pointer pKneeImage)
{
    mImage = pKneeImage;
}

bool CheckImageRod::Execute()
{
    ImageType::Pointer tImageClean = cleanImage(mImage, 0, 1000, -800, 10000);
    short threshold = 300;
    ImageType::Pointer otsuIm = OtsuMultipleThresholds(tImageClean, threshold);

    using BinaryThresholdFilterType = itk::BinaryThresholdImageFilter<ImageType, ImageType>;
    BinaryThresholdFilterType::Pointer BinaryFilter = BinaryThresholdFilterType::New();
    BinaryFilter->SetInput(otsuIm);

    const short  outsidevalue = 0;
    const short  insidevalue = 1000;

    BinaryFilter->SetOutsideValue(outsidevalue);
    BinaryFilter->SetInsideValue(insidevalue);

    const short lowerThreshold = -800;
    const short upperThreshold = 10000;

    BinaryFilter->SetUpperThreshold(upperThreshold);
    BinaryFilter->SetLowerThreshold(lowerThreshold);
    BinaryFilter->Update();// Running Filter;

    ////////////////////////////////////////////////////
    const unsigned int Dimension = 3;
    typedef unsigned short                                LabelType;
    typedef itk::Image< LabelType, Dimension >            OutputImageType;
    typedef itk::ShapeLabelObject< LabelType, Dimension > ShapeLabelObjectType;
    typedef itk::LabelMap< ShapeLabelObjectType >         LabelMapType;

    typedef itk::ConnectedComponentImageFilter <ImageType, OutputImageType >
        ConnectedComponentImageFilterType;
    typedef itk::LabelImageToShapeLabelMapFilter< OutputImageType, LabelMapType>
        I2LType;
    ConnectedComponentImageFilterType::Pointer connected =
        ConnectedComponentImageFilterType::New();
    connected->SetInput(BinaryFilter->GetOutput());
    connected->Update();
    typedef itk::LabelImageToShapeLabelMapFilter< OutputImageType, LabelMapType> I2LType;
    I2LType::Pointer i2l = I2LType::New();
    i2l->SetInput(connected->GetOutput());
    i2l->SetComputePerimeter(true);
    i2l->Update();
    LabelMapType *labelMap = i2l->GetOutput();

    LabelType rodLabel = 0;

    std::vector<LabelsInfo> LabelInfoVector;

    for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n)
    {
        ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);
        LabelInfoVector.push_back(LabelsInfo(labelObject->GetLabel(), labelObject->GetElongation(), labelObject->GetRegion()));
    }

    std::sort(LabelInfoVector.begin(), LabelInfoVector.end(), SortByElongation());

    /*for (int i = 0; i < 10; i++)
    {
        std::cout << LabelInfoVector[i].elongation <<std::endl;
    }*/

    ImageType::SizeType dims = otsuIm->GetLargestPossibleRegion().GetSize();

    if (LabelInfoVector.size() > 0)
    {
        if (LabelInfoVector[0].Similarity(dims) == false)
        {
            throw CheckImageRodException("The length of test rod does not have the size of the image.");
        }
        else
        {
            rodLabel = LabelInfoVector[0].MaxVolLabel;
        }
    }
    else
    {
        throw CheckImageRodException("It has not been possible to define a region for rod.");
    }

    LabelType* labelImageImBuffer = connected->GetOutput()->GetBufferPointer();
    ImageType::PixelType* myBuffer = otsuIm->GetBufferPointer();

    std::vector<itk::Point<double, 3>> centralPoints;
    double cont, x, y;

   /* ImageType::Pointer testImage = CloneImageTest<ImageType>(mImage);
    ImageType::PixelType* testBuffer = testImage->GetBufferPointer();*/

    for (int k = 0; k < dims[2]; k++)
    {
        cont = 0;
        x = 0;
        y = 0;

        for (int j = 0; j < dims[1]; j++)
        {
            for (int i = 0; i < dims[0]; i++)
            {
                int index = k * dims[1] * dims[0] + j * dims[0] + i;

                if (labelImageImBuffer[index] == rodLabel)
                {
                    cont++;
                    x = x + i;
                    y = y + j;
                }

                //if (labelImageImBuffer[index] == rodLabel)
                //{
                //    testBuffer[index] = 1;
                //}
                //else
                //{
                //    testBuffer[index] = 0;
                //}
            }
        }

        itk::Point<double, 3> tPoint;

        if (cont != 0)
        {
            x = x / cont;
            y = y / cont;
        }

        tPoint[0] = x;
        tPoint[1] = y;
        tPoint[2] = k;
        centralPoints.push_back(tPoint);
    }

    //SaveImageTest<ImageType>(testImage, "testRod");

    for (int i = 0; i < centralPoints.size() - 1; i++)
    {
        double diff = sqrt(pow((centralPoints[i][0] - centralPoints[i + 1][0]), 2) + pow((centralPoints[i][1] - centralPoints[i + 1][1]), 2) + pow((centralPoints[i][2] - centralPoints[i + 1][2]), 2));
        
        //std::cout << diff - 1. << std::endl;

        if (diff - 1. > 0.4)
        {
            throw CheckImageRodException("The rod found is not smooth.");
        }
    }

    return true;
}

ImageType::Pointer CheckImageRod::cleanImage(const ImageType::Pointer image, PixelType outsidevalue, PixelType insidevalue,
    PixelType lowerThreshold, PixelType upperThreshold)
{
    ImageType::Pointer result;

    using FilterType = itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType>;
    FilterType::Pointer smoothFilter = FilterType::New();

    smoothFilter->SetSigma(0.5);
    smoothFilter->SetInput(image);
    smoothFilter->Update();

    result = smoothFilter->GetOutput();

    using BinaryThresholdFilterType = itk::BinaryThresholdImageFilter<ImageType, ImageType>;
    BinaryThresholdFilterType::Pointer BinaryFilter = BinaryThresholdFilterType::New();
    BinaryFilter->SetInput(result);

    BinaryFilter->SetOutsideValue(outsidevalue);
    BinaryFilter->SetInsideValue(insidevalue);

    BinaryFilter->SetUpperThreshold(upperThreshold);
    BinaryFilter->SetLowerThreshold(lowerThreshold);
    BinaryFilter->Update();// Running Filter;

    ////////////////////////////////////////////////////
    const unsigned int Dimension = 3;
    typedef unsigned long long LabelType;
    typedef itk::Image< LabelType, Dimension >            OutputImageType;
    typedef itk::ShapeLabelObject< LabelType, Dimension > ShapeLabelObjectType;
    typedef itk::LabelMap< ShapeLabelObjectType >         LabelMapType;

    typedef itk::ConnectedComponentImageFilter<ImageType, OutputImageType >
        ConnectedComponentImageFilterType;
    typedef itk::LabelImageToShapeLabelMapFilter< OutputImageType, LabelMapType>
        I2LType;
    ConnectedComponentImageFilterType::Pointer connected =
        ConnectedComponentImageFilterType::New();

    connected->SetInput(BinaryFilter->GetOutput());
    connected->Update();
    typedef itk::LabelImageToShapeLabelMapFilter< OutputImageType, LabelMapType> I2LType;
    I2LType::Pointer i2l = I2LType::New();
    i2l->SetInput(connected->GetOutput());
    i2l->SetComputePerimeter(true);
    i2l->Update();
    LabelMapType *labelMap = i2l->GetOutput();

    LabelType maxVolLabel = 0;
    itk::SizeValueType maxVol = 0;

    for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n)
    {
        ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);

        if (labelObject->GetNumberOfPixels() > maxVol)
        {
            maxVolLabel = labelObject->GetLabel();
            maxVol = labelObject->GetNumberOfPixels();
        }
    }

    ImageType::RegionType inputRegion = result->GetLargestPossibleRegion();
    ImageType::SizeType dims = inputRegion.GetSize();
    LabelType* labelImageImBuffer = connected->GetOutput()->GetBufferPointer();

    ImageType::PixelType* imBuffer = result->GetBufferPointer();

    for (int k = 0; k < dims[2]; k++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int i = 0; i < dims[0]; i++)
            {
                int index = k * dims[1] * dims[0] + j * dims[0] + i;
                if (labelImageImBuffer[index] != maxVolLabel)
                {
                    imBuffer[index] = -1000;
                }
            }
        }
    }

    return result;
}

ImageType::Pointer CheckImageRod::OtsuMultipleThresholds(const ImageType::Pointer image, short &threshold)
{
    constexpr unsigned int Dimension = 3;
    using InputPixelType = short;
    using InternalPixelType = unsigned short;

    using InputImageType = itk::Image<InputPixelType, Dimension>;
    using InternalImageType = itk::Image<InternalPixelType, Dimension>;

    //using ReaderType = itk::ImageFileReader<InputImageType>;
    using OtsuMultipleThresholdsFilterType = itk::OtsuMultipleThresholdsImageFilter<ImageType, InternalImageType>;

    OtsuMultipleThresholdsFilterType::Pointer OtsuMultipleThresholdsfilter = OtsuMultipleThresholdsFilterType::New();
    //itk::SimpleFilterWatcher watcher(OtsuMultipleThresholdsfilter);

    auto numberOfHistogramBins = 100;
    OtsuMultipleThresholdsfilter->SetNumberOfHistogramBins(numberOfHistogramBins);
    //ITK_TEST_SET_GET_VALUE(numberOfHistogramBins, filter->GetNumberOfHistogramBins());

    auto numberOfThresholds = 4;
    OtsuMultipleThresholdsfilter->SetNumberOfThresholds(numberOfThresholds);
    //ITK_TEST_SET_GET_VALUE(numberOfThresholds, filter->GetNumberOfThresholds());

    auto labelOffset = 100;
    OtsuMultipleThresholdsfilter->SetLabelOffset(labelOffset);
    OtsuMultipleThresholdsfilter->SetInput(image);

    //ITK_TRY_EXPECT_NO_EXCEPTION(filter->Update());
    OtsuMultipleThresholdsfilter->Update();

    OtsuMultipleThresholdsFilterType::ThresholdVectorType thresholds = OtsuMultipleThresholdsfilter->GetThresholds();

    bool isHaveThreshold = false;
    for (unsigned int i = 0; i < thresholds.size(); i++)
    {
        //std::cout << thresholds[i] << std::endl;
        if (thresholds[i] > 200 && thresholds[i] < 500)
        {
            threshold = thresholds[i];
            isHaveThreshold = true;
            break;
        }
    }

    using BinaryThresholdFilterType = itk::BinaryThresholdImageFilter<ImageType, InternalImageType>;
    BinaryThresholdFilterType::Pointer BinaryFilter = BinaryThresholdFilterType::New();
    BinaryFilter->SetInput(image);

    const InternalPixelType  outsidevalue = -1000;
    const InternalPixelType  insidevalue = 1000;

    BinaryFilter->SetOutsideValue(outsidevalue);
    BinaryFilter->SetInsideValue(insidevalue);

    const InternalPixelType lowerThreshold = threshold;
    const InternalPixelType upperThreshold = 10000;

    BinaryFilter->SetUpperThreshold(upperThreshold);
    BinaryFilter->SetLowerThreshold(lowerThreshold);
    BinaryFilter->Update();// Running Filter;

    typedef itk::CastImageFilter<InternalImageType, ImageType > CastFilterType;
    CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput(BinaryFilter->GetOutput());
    castFilter->Update();
    return castFilter->GetOutput();
}