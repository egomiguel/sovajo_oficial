//#include "itkGDCMImageIO.h"
//#include "itkGDCMSeriesFileNames.h"
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

#include "vtkSimplePointsReader.h"
#include "vtkImageViewer2.h"
#include "vtkTextProperty.h"
#include "vtkTextMapper.h"

#include "vtkAutoInit.h"
#include "vtkSmartPointer.h"
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkCellPicker.h"
#include "vtkCommand.h"
#include "vtkImageActor.h"
#include "vtkImageReslice.h"
#include "vtkInteractorStyleImage.h"
#include "vtkImageMapToColors.h"
#include "vtkImagePlaneWidget.h"
#include "vtkImageReader.h"
#include "vtkInteractorEventRecorder.h"
#include "vtkLookupTable.h"
#include "vtkOutlineFilter.h"
#include "vtkDICOMImageReader.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkPlaneSource.h"
#include "vtkPlane.h"
#include "vtkResliceCursorActor.h"
#include "vtkResliceCursorPolyDataAlgorithm.h"
#include "vtkResliceCursor.h"
#include "vtkResliceCursorWidget.h"
#include "vtkResliceCursorLineRepresentation.h"
#include "vtkBiDimensionalWidget.h"

#include "vtkAxesActor.h"
#include "vtkTransform.h"
#include "vtkTextActor.h"
#include "vtkProperty2D.h"

#include "vtkImageToPolyDataFilter.h"
#include "itkOtsuMultipleThresholdsImageFilter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkCastImageFilter.h"
#include "itkExtractImageFilter.h"
#include "vtkClipPolyData.h"
#include "vtkCylinder.h"
#include "vtkImplicitBoolean.h"
#include "vtkImageDataGeometryFilter.h"

#include "vtkPolyData.h"
#include "vtkParametricSpline.h"
#include "vtkParametricFunctionSource.h"
#include "vtkSphereSource.h"
#include "vtkGlyph3DMapper.h"

#include "vtkCellArray.h"
#include "vtkInteractorStyleTrackballActor.h"
#include "vtkObjectFactory.h"

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
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include "SegmentationException.hpp"
#include "AutomaticSegmentation.hpp"
#include <fstream>
#include <stdio.h>

using namespace TKA::SEGMENTATION;

struct SLine
{
    cv::Point3d directVector, mPoint;

    SLine(const cv::Point3i& pA, const cv::Point3i& pB)
    {
        mPoint = cv::Point3d(double(pA.x), double(pA.y), double(pA.z));
        cv::Point3i temp = pA - pB;
        directVector = cv::Point3d(double(temp.x), double(temp.y), double(temp.z));
        //directVector = directVector / sqrt(directVector.dot(directVector));
    }

    cv::Point3d getProjectPoint(const cv::Point3d& pPoint) const
    {
        cv::Point3d diff = pPoint - mPoint;
        cv::Point3d projectOnDirector = ((diff.dot(directVector)) / (directVector.dot(directVector))) * directVector;
        cv::Point3d projection = mPoint + projectOnDirector;
        return projection;
    }
};

struct BoneReference
{
public:
    BoneReference()
    {
        femurShort = SegmentImageType::New();
        tibiaShort = SegmentImageType::New();
        kneeCapShort = SegmentImageType::New();
        fibulaShort = SegmentImageType::New();
    }

    void updateFemur(ImageType::IndexType index, ImageType::SizeType pSize)
    {
        SegmentImageType::IndexType start;
        start[0] = 0;  // first index on X
        start[1] = 0;  // first index on Y
        start[2] = 0;  // first index on Z

        /*std::vector<unsigned long long> v = { pSize[0], pSize[1], pSize[2] };

        unsigned long long min = *(std::min_element(v.begin(), v.end()));*/
        int min = 60;

        femurSize[0] = pSize[0] + min;
        femurSize[1] = pSize[1] + min;
        femurSize[2] = pSize[2] + min;

        SegmentImageType::RegionType region;
        region.SetSize(femurSize);
        region.SetIndex(start);

        femurShort->SetRegions(region);
        femurShort->Allocate();

        cv::Point3i intShort = cv::Point3i(min / 2, min / 2, min / 2);
        femurDiff = cv::Point3i(index[0], index[1], index[2]) - intShort;
    }

    void updateTibia(ImageType::IndexType index, ImageType::SizeType pSize)
    {
        SegmentImageType::IndexType start;
        start[0] = 0;  // first index on X
        start[1] = 0;  // first index on Y
        start[2] = 0;  // first index on Z

        /*std::vector<unsigned long long> v = { pSize[0], pSize[1], pSize[2] };
        unsigned long long min = *(std::min_element(v.begin(), v.end()));*/

        int min = 60;

        tibiaSize[0] = pSize[0] + min;
        tibiaSize[1] = pSize[1] + min;
        tibiaSize[2] = pSize[2] + min;

        SegmentImageType::RegionType region;
        region.SetSize(tibiaSize);
        region.SetIndex(start);

        tibiaShort->SetRegions(region);
        tibiaShort->Allocate();

        cv::Point3i intShort = cv::Point3i(min / 2, min / 2, min / 2);
        tibiaDiff = cv::Point3i(index[0], index[1], index[2]) - intShort;
    }

    void updateKnee(ImageType::IndexType index, ImageType::SizeType pSize)
    {
        SegmentImageType::IndexType start;
        start[0] = 0;  // first index on X
        start[1] = 0;  // first index on Y
        start[2] = 0;  // first index on Z

        int min = 60;

        kneeSize[0] = pSize[0] + min;
        kneeSize[1] = pSize[1] + min;
        kneeSize[2] = pSize[2] + min;

        SegmentImageType::RegionType region;
        region.SetSize(kneeSize);
        region.SetIndex(start);

        kneeCapShort->SetRegions(region);
        kneeCapShort->Allocate();

        cv::Point3i intShort = cv::Point3i(min / 2, min / 2, min / 2);
        kneeCapDiff = cv::Point3i(index[0], index[1], index[2]) - intShort;
    }

    void updateFibula(ImageType::IndexType index, ImageType::SizeType pSize)
    {
        SegmentImageType::IndexType start;
        start[0] = 0;  // first index on X
        start[1] = 0;  // first index on Y
        start[2] = 0;  // first index on Z

        int min = 60;

        fibulaSize[0] = pSize[0] + min;
        fibulaSize[1] = pSize[1] + min;
        fibulaSize[2] = pSize[2] + min;

        SegmentImageType::RegionType region;
        region.SetSize(fibulaSize);
        region.SetIndex(start);

        fibulaShort->SetRegions(region);
        fibulaShort->Allocate();

        cv::Point3i intShort = cv::Point3i(min / 2, min / 2, min / 2);
        fibulaDiff = cv::Point3i(index[0], index[1], index[2]) - intShort;
    }

    SegmentImageType::Pointer femurShort, tibiaShort, kneeCapShort, fibulaShort;
    SegmentImageType::SizeType femurSize, tibiaSize, kneeSize, fibulaSize;
    cv::Point3i femurDiff, tibiaDiff, kneeCapDiff, fibulaDiff;
};

struct LabelsInfo
{
    unsigned short MaxVolLabel;
    itk::SizeValueType pixels;
    double radius;
    ImageType::RegionType region;
    SegmentImageType::IndexType start;
    LabelsInfo(unsigned short MaxVolLabel, itk::SizeValueType pixels, ImageType::RegionType region)
    {
        this->MaxVolLabel = MaxVolLabel;
        this->pixels = pixels;
        this->region = region;
        ImageType::SizeType tSize = this->region.GetSize();
        start = this->region.GetIndex();
        if (tSize[0] >= tSize[1] && tSize[0] >= tSize[2])
        {
            this->radius = tSize[0];
        }

        else if (tSize[1] >= tSize[0] && tSize[1] >= tSize[2])
        {
            this->radius = tSize[1];
        }
        else
        {
            this->radius = tSize[2];
        }
    }

    bool Similarity(ImageType::SizeType pSize, bool strong = false)
    {
        ImageType::SizeType tSize = this->region.GetSize();

        int a = (pSize[0] - tSize[0]);
        int b = (pSize[1] - tSize[1]);
        int c = (pSize[2] - tSize[2]);

        int alpha = 10;

        if (strong == true)
        {
            alpha = 1;
        }

        if (abs(a) <= alpha || abs(b) <= alpha || abs(c) <= alpha)
        {
            return true;
        }
        return false;
    }
};

struct SortByPixels
{
    bool operator()(const LabelsInfo & label1, const LabelsInfo & label2) { return label1.pixels > label2.pixels; }
};

template<typename ImageType>
void SaveImageTest(typename ImageType::Pointer img, std::string name)
{
    std::cout << "Begin to save***" << std::endl;
    std::string fullName = name + ".nrrd";
    using WriterType = itk::ImageFileWriter<ImageType>;
    typename WriterType::Pointer writer = WriterType::New();

    writer->SetFileName(fullName);
    writer->SetInput(img);

    try
    {
        writer->Update();

    }
    catch (const itk::ExceptionObject & ex)
    {
        std::cout << ex.what() << std::endl;
    }
    std::cout << "Image saved" << std::endl;
}

VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);

AutomaticSegmentation::AutomaticSegmentation(ImageType::Pointer legImage)
{
    leg_image_ = legImage;
}

bool AutomaticSegmentation::ExecuteSegmentation()
{
    leg_image_clean_ = getOnlyLeg(leg_image_, 0, 1000, -800, 10000);

    //SaveImageTest<ImageType>(leg_image_clean_, "clean_image");

    short threshold = 300;
    ImageType::Pointer otsuIm = OtsuMultipleThresholds(leg_image_clean_, threshold); // Get threshold parameters

    //SaveImageTest<ImageType>(otsuIm, "otsu");

    getSegmentationBones(otsuIm);

    /*
    AutomaticSegmentation::BonesReference reference = getMarkersNew(otsuIm, leg_image_clean_);

    SaveImageTest<ImageType>(reference.marker, "markers");
    SaveImageTest<ImageType>(leg_image_clean_, "markers_cleanImage");

    femur_segment_ = WatershedSegmentation(leg_image_clean_, reference.marker);
    tibia_segment_ = CloneImage<SegmentImageType>(femur_segment_); // Clone the full image (complete leg)
    knee_cap_ = CloneImage<SegmentImageType>(femur_segment_);

    getBoneParts(femur_segment_, tibia_segment_, knee_cap_); // separate femur and tibia
    */

    /*SegmentImageType::Pointer tibiaIm = castImage(tibiaFirst);
    SegmentImageType::Pointer femoralIm = castImage(otsuIm);
    SegmentImageType::Pointer kneeCapIm = castImage(kneeCapFirst);

    femoralIm = BinaryDilateImage3d(femoralIm, 15);
    tibiaIm = BinaryDilateImage3d(tibiaIm, 20);
    kneeCapIm = BinaryDilateImage3d(kneeCapIm, 15);

    using FillholeFilterType = itk::BinaryFillholeImageFilter<SegmentImageType>;
    FillholeFilterType::Pointer fillHoleFemur = FillholeFilterType::New();
    FillholeFilterType::Pointer fillHoleTibia = FillholeFilterType::New();
    FillholeFilterType::Pointer fillHoleKneeCap = FillholeFilterType::New();

    fillHoleFemur->SetInput(femoralIm);
    fillHoleFemur->SetForegroundValue(1);
    fillHoleFemur->Update();

    fillHoleTibia->SetInput(tibiaIm);
    fillHoleTibia->SetForegroundValue(1);
    fillHoleTibia->Update();

    fillHoleKneeCap->SetInput(kneeCapIm);
    fillHoleKneeCap->SetForegroundValue(1);
    fillHoleKneeCap->Update();

    femur_segment_ = BinaryErodeImage3d(fillHoleFemur->GetOutput(), 15);
    tibia_segment_ = BinaryErodeImage3d(fillHoleTibia->GetOutput(), 20);
    knee_cap_ = BinaryErodeImage3d(fillHoleKneeCap->GetOutput(), 15);*/

    return true;
}

SegmentImageType::Pointer AutomaticSegmentation::closeBone(const SegmentImageType::Pointer image, short r) const
{
    //SegmentImageType::Pointer tIm = castImage(tibiaFirst);

    SegmentImageType::Pointer imageIn = BinaryDilateImage3d(image, r);

    using FillholeFilterType = itk::BinaryFillholeImageFilter<SegmentImageType>;
    FillholeFilterType::Pointer fillHole = FillholeFilterType::New();

    fillHole->SetInput(imageIn);
    fillHole->SetForegroundValue(1);
    fillHole->Update();

    return BinaryErodeImage3d(fillHole->GetOutput(), r);
}

SegmentImageType::Pointer AutomaticSegmentation::GetFemoralSegment() const
{
    return femur_segment_;
}

SegmentImageType::Pointer AutomaticSegmentation::GetTibiaSegment() const
{
    return tibia_segment_;
}

SegmentImageType::Pointer AutomaticSegmentation::GetKneeCapSegment() const
{
    return knee_cap_;
}

SegmentImageType::Pointer AutomaticSegmentation::GetFibulaSegment() const
{
    return fibula_segment;
}

SegmentImageType::Pointer AutomaticSegmentation::WatershedSegmentation(ImageType::Pointer image, ImageType::Pointer marker) const
{
    /*using FilterType = itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType>;
    FilterType::Pointer smoothFilter = FilterType::New();

    smoothFilter->SetSigma(0.5);
    smoothFilter->SetInput(image);
    smoothFilter->Update();*/

    //ImageType::Pointer result = image; // smoothFilter->GetOutput();

    //SaveImageTest<ImageType>(result, "image_today_gauss");

    typedef itk::BinaryBallStructuringElement< ImageType::PixelType, 3  > StructuringElementType;
    StructuringElementType  kernel;
    kernel.SetRadius(2);
    kernel.CreateStructuringElement();

    typedef itk::WhiteTopHatImageFilter< ImageType, ImageType, StructuringElementType > HatFilterType;

    HatFilterType::Pointer hatFilter = HatFilterType::New();
    hatFilter->SetInput(image);
    hatFilter->SetKernel(kernel);
    hatFilter->Update();

    ImageType::Pointer resultHat = hatFilter->GetOutput();

    SaveImageTest<ImageType>(resultHat, "image_today_hat");

    using OutputImageType = ImageType;
    using ConnectedComponentImageFilterType = itk::ConnectedComponentImageFilter<ImageType, OutputImageType>;
    ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
    connected->SetInput(marker);
    connected->SetFullyConnected(false);
    connected->Update();

    OutputImageType::Pointer connectedMarker = connected->GetOutput();

    using MorphologicalWatershedFilterType = itk::MorphologicalWatershedFromMarkersImageFilter<ImageType, OutputImageType>;
    MorphologicalWatershedFilterType::Pointer watershedFilter = MorphologicalWatershedFilterType::New();
    watershedFilter->SetInput1(resultHat);
    watershedFilter->SetInput2(connectedMarker);
    watershedFilter->SetFullyConnected(false);
    watershedFilter->SetMarkWatershedLine(true);
    watershedFilter->Update();
    OutputImageType::Pointer resultWater = watershedFilter->GetOutput();
    SaveImageTest<ImageType>(resultWater, "resultWater");

    using LabelType = PixelType;
    using LabelObjectType = itk::LabelObject<LabelType, 3>;
    using LabelMapType = itk::LabelMap<LabelObjectType>;

    using ConverterType = itk::LabelImageToLabelMapFilter<ImageType, LabelMapType>;
    ConverterType::Pointer converter = ConverterType::New();
    converter->SetInput(resultWater);
    converter->Update();

    LabelMapType *labelMap = converter->GetOutput();
    unsigned int tSize = labelMap->GetNumberOfLabelObjects();
    std::cout << "Number water: " << tSize << std::endl;
    PixelType bigLabel = 1;
    itk::SizeValueType regionSize = 0;

    for (unsigned int i = 0; i < tSize; ++i)
    {
        LabelMapType::LabelObjectType * labelObject = labelMap->GetNthLabelObject(i);
        if (labelObject->Size() > regionSize)
        {
            regionSize = labelObject->Size();
            bigLabel = labelObject->GetLabel();
        }
    }

    using ConnectedComponentImageFilterType2 = itk::ConnectedComponentImageFilter<ImageType, ImageType >;
    ConnectedComponentImageFilterType2::Pointer connected2 = ConnectedComponentImageFilterType2::New();
    connected2->SetInput(resultWater);
    connected2->Update();

    LabelType* labelImageImBuffer = connected2->GetOutput()->GetBufferPointer();

    ImageType::PixelType* fullBuffer = resultWater->GetBufferPointer();
    ImageType::SizeType dims = resultWater->GetLargestPossibleRegion().GetSize();

    for (int k = 0; k < dims[2]; k++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int i = 0; i < dims[0]; i++)
            {
                int index = k * dims[1] * dims[0] + j * dims[0] + i;
                if (labelImageImBuffer[index] == bigLabel)
                {
                    fullBuffer[index] = 0;
                }
                else
                {
                    fullBuffer[index] = 1;
                }
            }
        }
    }

    return castImage(resultWater);
}

void AutomaticSegmentation::getSegmentationBones(const ImageType::Pointer image)
{
    using BinaryThresholdFilterType = itk::BinaryThresholdImageFilter<ImageType, ImageType>;
    BinaryThresholdFilterType::Pointer BinaryFilter = BinaryThresholdFilterType::New();
    BinaryFilter->SetInput(image);

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

    LabelType femurLabel = 0;
    LabelType tibiaLabel = 0;
    LabelType kneeLabel = 0;
    LabelType fibulaLabel = 0;

    std::vector<LabelsInfo> LabelInfoVector, LabelInfoVectorTemp;

    for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n)
    {
        ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);
        LabelInfoVectorTemp.push_back(LabelsInfo(labelObject->GetLabel(), labelObject->GetNumberOfPixels(), labelObject->GetRegion()));
    }

    std::sort(LabelInfoVectorTemp.begin(), LabelInfoVectorTemp.end(), SortByPixels());

    ImageType::RegionType kneeRegion, fibulaRegion;

    if (LabelInfoVectorTemp.size() > 6)
    {
        ImageType::SizeType dims = image->GetLargestPossibleRegion().GetSize();

        for (int i = 0; i < 7; i++)
        {
            if (LabelInfoVectorTemp[i].Similarity(dims) == false)
            {
                LabelInfoVector.push_back(LabelInfoVectorTemp[i]);
            }
        }

        if (LabelInfoVector.size() >= 4)
        {
            if (LabelInfoVector[0].start[2] > LabelInfoVector[1].start[2])
            {
                femurLabel = LabelInfoVector[0].MaxVolLabel;
                tibiaLabel = LabelInfoVector[1].MaxVolLabel;
            }
            else
            {
                femurLabel = LabelInfoVector[1].MaxVolLabel;
                tibiaLabel = LabelInfoVector[0].MaxVolLabel;
            }

            if (LabelInfoVector[2].radius < LabelInfoVector[3].radius)
            {
                kneeLabel = LabelInfoVector[2].MaxVolLabel;
                kneeRegion = LabelInfoVector[2].region;

                fibulaLabel = LabelInfoVector[3].MaxVolLabel;
                fibulaRegion = LabelInfoVector[3].region;
            }
            else
            {
                kneeLabel = LabelInfoVector[3].MaxVolLabel;
                kneeRegion = LabelInfoVector[3].region;

                fibulaLabel = LabelInfoVector[2].MaxVolLabel;
                fibulaRegion = LabelInfoVector[2].region;
            }
        }
        else
        {
            throw SegmentationExceptionCode::CAN_NOT_DEFINE_SEGMENTATION_REGION_FOR_EACH_BONE;
        }
    }
    else
    {
        throw SegmentationExceptionCode::CAN_NOT_DEFINE_SEGMENTATION_REGION_FOR_EACH_BONE;
    }

    ImageType::RegionType femurRegion = LabelInfoVector[0].region;
    ImageType::RegionType tibiaRegion = LabelInfoVector[1].region;

    ImageType::IndexType tibiaCorner = tibiaRegion.GetIndex();
    ImageType::SizeType tibiaSize = tibiaRegion.GetSize();

    ImageType::IndexType kneeCapCorner = kneeRegion.GetIndex();
    ImageType::SizeType kneeCapSize = kneeRegion.GetSize();

    ImageType::IndexType femurCorner = femurRegion.GetIndex();
    ImageType::SizeType femurSize = femurRegion.GetSize();

    ImageType::IndexType fibulaCorner = fibulaRegion.GetIndex();
    ImageType::SizeType fibulaSize = fibulaRegion.GetSize();

    ImageType::RegionType inputRegion = image->GetLargestPossibleRegion();
    ImageType::SizeType dims = inputRegion.GetSize();
    LabelType* labelImageImBuffer = connected->GetOutput()->GetBufferPointer();

    BoneReference reference;
    reference.updateFemur(femurCorner, femurSize);
    reference.updateTibia(tibiaCorner, tibiaSize);
    reference.updateKnee(kneeCapCorner, kneeCapSize);
    reference.updateFibula(fibulaCorner, fibulaSize);

    SegmentImageType::PixelType* kneeBuffer = reference.kneeCapShort->GetBufferPointer();
    SegmentImageType::PixelType* femurBuffer = reference.femurShort->GetBufferPointer();
    SegmentImageType::PixelType* tibiaBuffer = reference.tibiaShort->GetBufferPointer();
    SegmentImageType::PixelType* fibulaBuffer = reference.fibulaShort->GetBufferPointer();

    for (int k = 0; k < reference.kneeSize[2]; k++)
    {
        for (int j = 0; j < reference.kneeSize[1]; j++)
        {
            for (int i = 0; i < reference.kneeSize[0]; i++)
            {
                cv::Point3i mainCoord = reference.kneeCapDiff + cv::Point3i(i, j, k);
                int indexShort = k * reference.kneeSize[1] * reference.kneeSize[0] + j * reference.kneeSize[0] + i;
                int indexBig = mainCoord.z * dims[1] * dims[0] + mainCoord.y * dims[0] + mainCoord.x;

                if (mainCoord.x >= dims[0] || mainCoord.y >= dims[1] || mainCoord.z >= dims[2])
                {
                    kneeBuffer[indexShort] = 0;
                    continue;
                }

                if (labelImageImBuffer[indexBig] != kneeLabel)
                {
                    kneeBuffer[indexShort] = 0;
                }
                else
                {
                    kneeBuffer[indexShort] = 1;
                }
            }
        }
    }

    for (int k = 0; k < reference.femurSize[2]; k++)
    {
        for (int j = 0; j < reference.femurSize[1]; j++)
        {
            for (int i = 0; i < reference.femurSize[0]; i++)
            {
                cv::Point3i mainCoord = reference.femurDiff + cv::Point3i(i, j, k);
                int indexShort = k * reference.femurSize[1] * reference.femurSize[0] + j * reference.femurSize[0] + i;
                int indexBig = mainCoord.z * dims[1] * dims[0] + mainCoord.y * dims[0] + mainCoord.x;

                if (mainCoord.x >= dims[0] || mainCoord.y >= dims[1] || mainCoord.z >= dims[2])
                {
                    femurBuffer[indexShort] = 0;
                    continue;
                }

                if (labelImageImBuffer[indexBig] == femurLabel)
                {
                    femurBuffer[indexShort] = 1;
                }
                else
                {
                    femurBuffer[indexShort] = 0;
                }
            }
        }
    }

    for (int k = 0; k < reference.tibiaSize[2]; k++)
    {
        for (int j = 0; j < reference.tibiaSize[1]; j++)
        {
            for (int i = 0; i < reference.tibiaSize[0]; i++)
            {
                cv::Point3i mainCoord = reference.tibiaDiff + cv::Point3i(i, j, k);
                int indexShort = k * reference.tibiaSize[1] * reference.tibiaSize[0] + j * reference.tibiaSize[0] + i;
                int indexBig = mainCoord.z * dims[1] * dims[0] + mainCoord.y * dims[0] + mainCoord.x;

                if (mainCoord.x >= dims[0] || mainCoord.y >= dims[1] || mainCoord.z >= dims[2])
                {
                    tibiaBuffer[indexShort] = 0;
                    continue;
                }

                if (labelImageImBuffer[indexBig] == tibiaLabel)
                {
                    tibiaBuffer[indexShort] = 1;
                }
                else
                {
                    tibiaBuffer[indexShort] = 0;
                }
            }
        }
    }

    for (int k = 0; k < reference.fibulaSize[2]; k++)
    {
        for (int j = 0; j < reference.fibulaSize[1]; j++)
        {
            for (int i = 0; i < reference.fibulaSize[0]; i++)
            {
                cv::Point3i mainCoord = reference.fibulaDiff + cv::Point3i(i, j, k);
                int indexShort = k * reference.fibulaSize[1] * reference.fibulaSize[0] + j * reference.fibulaSize[0] + i;
                int indexBig = mainCoord.z * dims[1] * dims[0] + mainCoord.y * dims[0] + mainCoord.x;

                if (mainCoord.x >= dims[0] || mainCoord.y >= dims[1] || mainCoord.z >= dims[2])
                {
                    fibulaBuffer[indexShort] = 0;
                    continue;
                }

                if (labelImageImBuffer[indexBig] == fibulaLabel)
                {
                    fibulaBuffer[indexShort] = 1;
                }
                else
                {
                    fibulaBuffer[indexShort] = 0;
                }
            }
        }
    }

    SegmentImageType::Pointer femurClose = closeBone(reference.femurShort, 20);
    SegmentImageType::Pointer tibiaClose = closeBone(reference.tibiaShort, 20);
    SegmentImageType::Pointer kneeClose = closeBone(reference.kneeCapShort, 15);
    SegmentImageType::Pointer fibulaClose = closeBone(reference.fibulaShort, 15);

    femur_segment_ = castImage(image);
    tibia_segment_ = castImage(image);
    knee_cap_ = castImage(image);
    fibula_segment = castImage(image);

    SegmentImageType::RegionType myFemurRegion = femur_segment_->GetLargestPossibleRegion();
    SegmentImageType::SizeType myFemurDims = myFemurRegion.GetSize();
    SegmentImageType::PixelType* myFemurBuffer = femur_segment_->GetBufferPointer();
    SegmentImageType::PixelType* myFemurBufferShort = femurClose->GetBufferPointer();

    for (int k = 0; k < myFemurDims[2]; k++)
    {
        for (int j = 0; j < myFemurDims[1]; j++)
        {
            for (int i = 0; i < myFemurDims[0]; i++)
            {
                int indexBig = k * myFemurDims[1] * myFemurDims[0] + j * myFemurDims[0] + i;
                cv::Point3i mainCoord = cv::Point3i(i, j, k) - reference.femurDiff;

                if (mainCoord.x >= reference.femurSize[0] || mainCoord.y >= reference.femurSize[1] || mainCoord.z >= reference.femurSize[2] ||
                    mainCoord.x < 0 || mainCoord.y < 0 || mainCoord.z < 0)
                {
                    myFemurBuffer[indexBig] = 0;
                    continue;
                }

                int indexShort = mainCoord.z * reference.femurSize[1] * reference.femurSize[0] + mainCoord.y * reference.femurSize[0] + mainCoord.x;

                myFemurBuffer[indexBig] = myFemurBufferShort[indexShort];
            }
        }
    }

    SegmentImageType::RegionType myTibiaRegion = tibia_segment_->GetLargestPossibleRegion();
    SegmentImageType::SizeType myTibiaDims = myTibiaRegion.GetSize();
    SegmentImageType::PixelType* myTibiaBuffer = tibia_segment_->GetBufferPointer();
    SegmentImageType::PixelType* myTibiaBufferShort = tibiaClose->GetBufferPointer();

    for (int k = 0; k < myTibiaDims[2]; k++)
    {
        for (int j = 0; j < myTibiaDims[1]; j++)
        {
            for (int i = 0; i < myTibiaDims[0]; i++)
            {
                int indexBig = k * myTibiaDims[1] * myTibiaDims[0] + j * myTibiaDims[0] + i;
                cv::Point3i mainCoord = cv::Point3i(i, j, k) - reference.tibiaDiff;

                if (mainCoord.x >= reference.tibiaSize[0] || mainCoord.y >= reference.tibiaSize[1] || mainCoord.z >= reference.tibiaSize[2] ||
                    mainCoord.x < 0 || mainCoord.y < 0 || mainCoord.z < 0)
                {
                    myTibiaBuffer[indexBig] = 0;
                    continue;
                }

                int indexShort = mainCoord.z * reference.tibiaSize[1] * reference.tibiaSize[0] + mainCoord.y * reference.tibiaSize[0] + mainCoord.x;

                myTibiaBuffer[indexBig] = myTibiaBufferShort[indexShort];
            }
        }
    }

    SegmentImageType::RegionType myKneeRegion = knee_cap_->GetLargestPossibleRegion();
    SegmentImageType::SizeType myKneeDims = myKneeRegion.GetSize();
    SegmentImageType::PixelType* myKneeBuffer = knee_cap_->GetBufferPointer();
    SegmentImageType::PixelType* myKneeBufferShort = kneeClose->GetBufferPointer();

    for (int k = 0; k < myKneeDims[2]; k++)
    {
        for (int j = 0; j < myKneeDims[1]; j++)
        {
            for (int i = 0; i < myKneeDims[0]; i++)
            {
                int indexBig = k * myKneeDims[1] * myKneeDims[0] + j * myKneeDims[0] + i;
                cv::Point3i mainCoord = cv::Point3i(i, j, k) - reference.kneeCapDiff;

                if (mainCoord.x >= reference.kneeSize[0] || mainCoord.y >= reference.kneeSize[1] || mainCoord.z >= reference.kneeSize[2] ||
                    mainCoord.x < 0 || mainCoord.y < 0 || mainCoord.z < 0)
                {
                    myKneeBuffer[indexBig] = 0;
                    continue;
                }

                int indexShort = mainCoord.z * reference.kneeSize[1] * reference.kneeSize[0] + mainCoord.y * reference.kneeSize[0] + mainCoord.x;

                myKneeBuffer[indexBig] = myKneeBufferShort[indexShort];
            }
        }
    }

    SegmentImageType::RegionType myFibulaRegion = fibula_segment->GetLargestPossibleRegion();
    SegmentImageType::SizeType myFibulaDims = myFibulaRegion.GetSize();
    SegmentImageType::PixelType* myFibulaBuffer = fibula_segment->GetBufferPointer();
    SegmentImageType::PixelType* myFibulaBufferShort = fibulaClose->GetBufferPointer();

    for (int k = 0; k < myFibulaDims[2]; k++)
    {
        for (int j = 0; j < myFibulaDims[1]; j++)
        {
            for (int i = 0; i < myFibulaDims[0]; i++)
            {
                int indexBig = k * myFibulaDims[1] * myFibulaDims[0] + j * myFibulaDims[0] + i;
                cv::Point3i mainCoord = cv::Point3i(i, j, k) - reference.fibulaDiff;

                if (mainCoord.x >= reference.fibulaSize[0] || mainCoord.y >= reference.fibulaSize[1] || mainCoord.z >= reference.fibulaSize[2] ||
                    mainCoord.x < 0 || mainCoord.y < 0 || mainCoord.z < 0)
                {
                    myFibulaBuffer[indexBig] = 0;
                    continue;
                }

                int indexShort = mainCoord.z * reference.fibulaSize[1] * reference.fibulaSize[0] + mainCoord.y * reference.fibulaSize[0] + mainCoord.x;

                myFibulaBuffer[indexBig] = myFibulaBufferShort[indexShort];
            }
        }
    }
}

/*
ImageType::Pointer AutomaticSegmentation::getMarkers(const ImageType::Pointer image, ImageType::Pointer cleanImage) const
{
    using BinaryThresholdFilterType = itk::BinaryThresholdImageFilter<ImageType, ImageType>;
    BinaryThresholdFilterType::Pointer BinaryFilter = BinaryThresholdFilterType::New();
    BinaryFilter->SetInput(image);

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

    LabelType femurLabel = 0;
    LabelType tibiaLabel = 0;
    LabelType kneeLabel = 0;

    std::vector<LabelsInfo> LabelInfoVector;

    for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n)
    {
        ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);
        LabelInfoVector.push_back(LabelsInfo(labelObject->GetLabel(), labelObject->GetNumberOfPixels(), labelObject->GetRoundness(), labelObject->GetElongation(), labelObject->GetRegion()));
    }

    std::sort(LabelInfoVector.begin(), LabelInfoVector.end(), SortByPixels());

    ImageType::RegionType kneeRegion;

    if (LabelInfoVector.size() >= 4)
    {
        femurLabel = LabelInfoVector[0].MaxVolLabel;
        tibiaLabel = LabelInfoVector[1].MaxVolLabel;
        if (LabelInfoVector[2].radius < LabelInfoVector[3].radius)
        {
            kneeLabel = LabelInfoVector[2].MaxVolLabel;
            kneeRegion = LabelInfoVector[2].region;
        }
        else
        {
            kneeLabel = LabelInfoVector[3].MaxVolLabel;
            kneeRegion = LabelInfoVector[3].region;
        }
    }
    else
    {
        throw SegmentationException("It has not been possible to define a region for each bone.");
    }

    ImageType::RegionType femurRegion = LabelInfoVector[0].region;
    ImageType::RegionType tibiaRegion = LabelInfoVector[1].region;

    ImageType::IndexType tibiaCorner = tibiaRegion.GetIndex();
    ImageType::SizeType tibiaSize = tibiaRegion.GetSize();

    ImageType::IndexType kneeCapCorner = kneeRegion.GetIndex();
    ImageType::SizeType kneeCapSize = kneeRegion.GetSize();

    ImageType::IndexType femurCorner = femurRegion.GetIndex();
    ImageType::SizeType femurSize = femurRegion.GetSize();

    cv::Point3i centerTibia = (cv::Point3i(tibiaCorner[0], tibiaCorner[1], tibiaCorner[2]) + cv::Point3i(tibiaSize[0], tibiaSize[1], tibiaSize[2])) / 2;

    cv::Point3i initKneeCap = cv::Point3i(kneeCapCorner[0], kneeCapCorner[1], kneeCapCorner[2]);

    cv::Point3i centerFemur = (cv::Point3i(femurCorner[0], femurCorner[1], femurCorner[2]) + cv::Point3i(femurSize[0], femurSize[1], femurSize[2])) / 2;

    int tibiaMaxZ = tibiaCorner[2] + tibiaSize[2] - (tibiaSize[2] / 5);
    int tibiaMinZ = tibiaCorner[2] + (tibiaSize[2] / 5);

    int femurMaxZ = femurCorner[2] + femurSize[2] - (femurSize[2] / 5);
    int femurMinZ = femurCorner[2] + (femurSize[2] / 4);

    int patellaMaxZ = kneeCapCorner[2] + kneeCapSize[2] - (kneeCapSize[2] / 4);
    int patellaMinZ = kneeCapCorner[2] + (kneeCapSize[2] / 4);

    int tibiaCont = 0;
    cv::Point3i tibiaPoint = cv::Point3i(0, 0, 0);

    int femurCont = 0;
    cv::Point3i femurPoint = cv::Point3i(0, 0, 0);

    int patellaCont = 0;
    cv::Point3i patellaPoint = cv::Point3i(0, 0, 0);

    ImageType::RegionType inputRegion = image->GetLargestPossibleRegion();
    ImageType::SizeType dims = inputRegion.GetSize();
    LabelType* labelImageImBuffer = connected->GetOutput()->GetBufferPointer();

    ImageType::Pointer fullImage = CloneImage<ImageType>(leg_image_clean_);
    ImageType::PixelType* fullBuffer = fullImage->GetBufferPointer();

    ImageType::PixelType* fullBufferClean = cleanImage->GetBufferPointer();

    for (int k = 0; k < dims[2]; k++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int i = 0; i < dims[0]; i++)
            {
                int index = k * dims[1] * dims[0] + j * dims[0] + i;

                ImageType::IndexType indexPos;
                indexPos[0] = i;
                indexPos[1] = j;
                indexPos[2] = k;

                if (femurRegion.IsInside(indexPos) == false && tibiaRegion.IsInside(indexPos) == false && kneeRegion.IsInside(indexPos) == false)
                {
                    fullBuffer[index] = 1;
                    fullBufferClean[index] = 0;
                }
                else
                {

                    fullBuffer[index] = 0;

                }

                if (labelImageImBuffer[index] == femurLabel)
                {
                    //fullBufferClean[index] = 1000;
                    if (femurMinZ == k)
                    {
                        femurPoint = femurPoint + cv::Point3i(i, j, k);
                        femurCont++;

                    }
                    else if (femurMinZ <= femurMaxZ && femurMinZ + 1 == k)
                    {
                        femurPoint = femurPoint / femurCont;

                        for (int z = -10; z <= 10; z++)
                        {
                            cv::Point3i myTemp = femurPoint + z * cv::Point3i(1, 1, 0);
                            int myIndex = myTemp.z * dims[1] * dims[0] + myTemp.y * dims[0] + myTemp.x;
                            fullBuffer[myIndex] = 2;
                        }

                        femurPoint.x = i;
                        femurPoint.y = j;
                        femurPoint.z = k;

                        femurCont = 1;

                        if (femurMinZ < femurMaxZ)
                        {
                            femurMinZ = k;
                        }
                        else
                        {
                            femurMinZ = -2;
                        }
                    }

                }

                else if (labelImageImBuffer[index] == tibiaLabel)
                {
                    //fullBufferClean[index] = 1000;

                    if (tibiaMinZ == k)
                    {
                        tibiaPoint = tibiaPoint + cv::Point3i(i, j, k);
                        tibiaCont++;

                    }
                    else if (tibiaMinZ <= tibiaMaxZ && tibiaMinZ + 1 == k)
                    {

                        tibiaPoint = tibiaPoint / tibiaCont;

                        for (int z = -10; z <= 10; z++)
                        {
                            cv::Point3i myTemp = tibiaPoint + z * cv::Point3i(1, 1, 0);
                            int myIndex = myTemp.z * dims[1] * dims[0] + myTemp.y * dims[0] + myTemp.x;
                            fullBuffer[myIndex] = 3;
                        }

                        tibiaPoint.x = i;
                        tibiaPoint.y = j;
                        tibiaPoint.z = k;

                        tibiaCont = 1;

                        if (tibiaMinZ < tibiaMaxZ)
                        {
                            tibiaMinZ = k;
                        }
                        else
                        {
                            tibiaMinZ = -2;
                        }
                    }
                }

                else if (labelImageImBuffer[index] == kneeLabel)
                {
                    //fullBufferClean[index] = 1000;

                    if (patellaMinZ == k)
                    {
                        patellaPoint = patellaPoint + cv::Point3i(i, j, k);
                        patellaCont++;

                    }
                    else if (patellaMinZ <= patellaMaxZ && patellaMinZ + 1 == k)
                    {

                        patellaPoint = patellaPoint / patellaCont;

                        for (int z = -5; z <= 5; z++)
                        {
                            cv::Point3i myTemp = patellaPoint + z * cv::Point3i(1, 1, 0);
                            int myIndex = myTemp.z * dims[1] * dims[0] + myTemp.y * dims[0] + myTemp.x;
                            fullBuffer[myIndex] = 4;
                        }

                        patellaPoint.x = i;
                        patellaPoint.y = j;
                        patellaPoint.z = k;

                        patellaCont = 1;

                        if (patellaMinZ < patellaMaxZ)
                        {
                            patellaMinZ = k;
                        }
                        else
                        {
                            patellaMinZ = -2;
                        }
                    }

                }
                //else if(!(femurRegion.IsInside(indexPos) == false && tibiaRegion.IsInside(indexPos) == false && kneeRegion.IsInside(indexPos) == false))
                //{
                //    fullBufferClean[index] = 500;
                //}
                else
                {
                    //fullBufferClean[index] = 0;
                }
            }
        }
    }

    SaveImageTest<ImageType>(fullImage, "image_today_markers");

    return fullImage;
}
*/

template<typename ImageType>
typename ImageType::Pointer AutomaticSegmentation::CloneImage(const typename ImageType::Pointer input) const
{
    using DuplicatorType = itk::ImageDuplicator<ImageType>;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(input);
    duplicator->Update();
    return duplicator->GetOutput();
}

/*
ImageType::Pointer AutomaticSegmentation::ExecuteSegmentation(const char* wFemoralDir, const char* wTibiaDir, const std::string& seriesName)
{
    using PixelType = signed short;

    using NamesGeneratorType = itk::GDCMSeriesFileNames;
    NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

    nameGenerator->SetUseSeriesDetails(true);
    nameGenerator->AddSeriesRestriction("0008|0021");
    nameGenerator->SetGlobalWarningDisplay(false);
    nameGenerator->SetDirectory(dirName);

    try {
        using SeriesIdContainer = std::vector<std::string>;
        const SeriesIdContainer& seriesUID = nameGenerator->GetSeriesUIDs();
        auto seriesItr = seriesUID.begin();
        auto seriesEnd = seriesUID.end();

        if (seriesItr != seriesEnd) {
            //std::cout << "The directory: ";
            //std::cout << dirName << std::endl;
            //std::cout << "Contains the following DICOM Series: ";
            //std::cout << std::endl;
        }
        else {
            //std::cout << "No DICOMs in: " << dirName << std::endl;
            return nullptr;
        }
        while (seriesItr != seriesEnd) {
            //std::cout << seriesItr->c_str() << std::endl;
            ++seriesItr;
        }
        seriesItr = seriesUID.begin();
        while (seriesItr != seriesUID.end()) {
            //std::cout << "---------------------------------------------" << std::endl;
            std::string seriesIdentifier;
            if (!seriesName.empty()) // If seriesIdentifier given convert only that
            {
                seriesIdentifier = seriesName;
                seriesItr = seriesUID.end();
            }
            else // otherwise convert everything
            {
                seriesIdentifier = seriesItr->c_str();
                seriesItr++;
            }
            //std::cout << "\nReading: ";
            //std::cout << seriesIdentifier << std::endl;
            using FileNamesContainer = std::vector<std::string>;
            FileNamesContainer fileNames = nameGenerator->GetFileNames(seriesIdentifier);
            if (fileNames.empty()) {
                continue;
            }

            using ReaderType = itk::ImageSeriesReader<ImageType>;
            ReaderType::Pointer reader = ReaderType::New();
            using ImageIOType = itk::GDCMImageIO;
            ImageIOType::Pointer dicomIO = ImageIOType::New();
            reader->SetImageIO(dicomIO);
            reader->SetFileNames(fileNames);
            // reader->ForceOrthogonalDirectionOff(); // properly read CTs with gantry tilt

            reader->Update();
            // Read only first series
            //typename ImageType::Pointer res = reader->GetOutput();
            ImageType::Pointer res = reader->GetOutput();
            /////////////////????//////////////////////////////
            using WriterType = itk::ImageFileWriter<ImageType>;
            WriterType::Pointer writer = WriterType::New();

            writer->SetFileName("leg.nrrd");
            writer->SetInput(res);

            std::cout << "Writing the image " << std::endl << std::endl;
            try
            {
                writer->Update();
            }
            catch (const itk::ExceptionObject & ex)
            {
                std::cout << ex.what() << std::endl;
            }
            /////////////////////////////////////////////////

            //ImageType::Pointer onlyLeg = getOnlyLeg<ImageType, ImageType::ImageType>(res, 0, 1000, -800, 10000);
            ImageType::Pointer onlyLeg = getOnlyLeg(res, 0, 1000, -800, 10000);
            short threshold = 300;
            //typename ImageType::Pointer otsuIm = OtsuMultipleThresholds<ImageType>(onlyLeg, threshold);
            ImageType::Pointer otsuIm = OtsuMultipleThresholds(onlyLeg, threshold);
            //vis3d(otsuIm);
            ImageType::Pointer tibiaIm = ImageType::New();
            //DeepCopy<ImageType>(otsuIm, tibiaIm);
            DeepCopy(otsuIm, tibiaIm);
            //typename ImageType::Pointer femoralIm = getTwoBone<ImageType>(otsuIm, tibiaIm);
            ImageType::Pointer femoralIm = getTwoBone(otsuIm, tibiaIm);
            //vis3d(femoralIm);
            //femoralIm = BinaryDilateImage3d<ImageType>(femoralIm, 25);
            //tibiaIm = BinaryDilateImage3d<ImageType>(tibiaIm, 25);

            femoralIm = BinaryDilateImage3d(femoralIm, 25);
            tibiaIm = BinaryDilateImage3d(tibiaIm, 25);

            femoralIm = BinaryErodeImage3d(femoralIm, 22);
            tibiaIm = BinaryErodeImage3d(tibiaIm, 22);
            //femoralIm = copyIm_byMark<ImageType>(res, femoralIm, 1000);
            //tibiaIm = copyIm_byMark<ImageType>(res, tibiaIm, 1000);

            SegmentImageType::Pointer FinalFemoral = copyIm_byMark(res, femoralIm, 1000);
            SegmentImageType::Pointer FinalTibia = copyIm_byMark(res, tibiaIm, 1000);

            //////////////////////////////////////////////////////////////////////
            using WriterType1 = itk::ImageFileWriter<SegmentImageType>;
            WriterType1::Pointer writer1 = WriterType1::New();

            writer1->SetFileName("segmentation3.nrrd");
            writer1->SetInput(FinalFemoral);

            std::cout << "Writing the image " << std::endl << std::endl;
            try
            {
                writer1->Update();
            }
            catch (const itk::ExceptionObject & ex)
            {
                std::cout << ex.what() << std::endl;
            }
            //////////////////////////////////////////////////////////////

            //femoralIm = copyIm_byMark(res, femoralIm, 1000);
            //tibiaIm = copyIm_byMark(res, tibiaIm, 1000);
            //////////////////////////////
            const char * outputDirectory = wFemoralDir;
            itksys::SystemTools::MakeDirectory(outputDirectory);
            using OutputPixelType = PixelType;
            constexpr unsigned int OutputDimension = 2;
            using Image2DType = itk::Image<OutputPixelType, OutputDimension>;
            using SeriesWriterType = itk::ImageSeriesWriter<ImageType, Image2DType>;
            SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
            seriesWriter->SetInput(femoralIm);
            seriesWriter->SetImageIO(dicomIO);
            nameGenerator->SetOutputDirectory(outputDirectory);
            NamesGeneratorType::FileNamesContainerType names = nameGenerator->GetOutputFileNames();

            ReaderType::DictionaryArrayRawPointer dtryArray = reader->GetMetaDataDictionaryArray();
            seriesWriter->SetFileNames(names);
            seriesWriter->SetMetaDataDictionaryArray(reader->GetMetaDataDictionaryArray());
            seriesWriter->Update();
            ////////////////////////?????////////////////////////
            outputDirectory = wTibiaDir;
            itksys::SystemTools::MakeDirectory(outputDirectory);
            SeriesWriterType::Pointer tibiaImSeriesWriter = SeriesWriterType::New();
            tibiaImSeriesWriter->SetInput(tibiaIm);
            tibiaImSeriesWriter->SetImageIO(dicomIO);
            nameGenerator->SetOutputDirectory(outputDirectory);
            NamesGeneratorType::FileNamesContainerType tibiaImNames = nameGenerator->GetOutputFileNames();


            tibiaImSeriesWriter->SetFileNames(tibiaImNames);
            tibiaImSeriesWriter->SetMetaDataDictionaryArray(reader->GetMetaDataDictionaryArray());
            tibiaImSeriesWriter->Update();
            std::cout << "----------0000-------------------11111-----------------------2222-----------------" << std::endl;
            return res;
        }
    }
    catch (itk::ExceptionObject& ex) {
        std::cout << ex << std::endl;
        return nullptr;

    }
    return nullptr;
}
*/

double AutomaticSegmentation::getSmoothParameter(const ImageType::Pointer image, PixelType outsidevalue, PixelType insidevalue,
    PixelType lowerThreshold, PixelType upperThreshold) const
{
    using BinaryThresholdFilterType = itk::BinaryThresholdImageFilter<ImageType, ImageType>;
    BinaryThresholdFilterType::Pointer BinaryFilter = BinaryThresholdFilterType::New();
    BinaryFilter->SetInput(image);

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

    unsigned long long tSize = labelMap->GetNumberOfLabelObjects();

    if (tSize >= 1000000)
    {
        return 0.8;
    }
    else if (tSize >= 850000)
    {
        return 0.7;
    }
    else if (tSize >= 700000)
    {
        return 0.6;
    }
    else if (tSize >= 550000)
    {
        return 0.5;
    }
    else if (tSize >= 400000)
    {
        return 0.4;
    }
    else if (tSize >= 250000)
    {
        return 0.3;
    }
    else if (tSize >= 150000)
    {
        return 0.25;
    }
    else if (tSize >= 70000)
    {
        return 0.2;
    }
    else if (tSize >= 36000)
    {
        return 0.1;
    }
    else
    {
        return 0;
    }
}

ImageType::Pointer AutomaticSegmentation::GradientAnisotropicDiffusion(const ImageType::Pointer image)
{
    const int numberOfIterations = 5;
    const float timeStep = 0.06;
    const float conductance = 2.0;

    constexpr unsigned int Dimension = 3;
    using OutPixelType = float;
    using OutputImageType = itk::Image<OutPixelType, Dimension>;

    using FilterCast1 = itk::CastImageFilter<ImageType, OutputImageType>;
    FilterCast1::Pointer filterCast1 = FilterCast1::New();

    filterCast1->SetInput(image);
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

ImageType::Pointer AutomaticSegmentation::getOnlyLeg(const ImageType::Pointer image, PixelType outsidevalue, PixelType insidevalue,
    PixelType lowerThreshold, PixelType upperThreshold)
{
    ImageType::Pointer result;
    
    /*
    double tParameter = getSmoothParameter(image, outsidevalue, insidevalue, lowerThreshold, upperThreshold);
    std::cout << "Filter: " << tParameter << std::endl;
    if (tParameter > 0.1)
    {
        using FilterType = itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType>;
        FilterType::Pointer smoothFilter = FilterType::New();

        smoothFilter->SetSigma(tParameter);
        smoothFilter->SetInput(image);
        smoothFilter->Update();

        result = smoothFilter->GetOutput();
    }
    else
    {
        result = image;
    }
    */
    result = GradientAnisotropicDiffusion(image);

    //SaveImageTest<ImageType>(result, "denoising");

    using BinaryThresholdFilterType = itk::BinaryThresholdImageFilter<ImageType, ImageType>;
    BinaryThresholdFilterType::Pointer BinaryFilter = BinaryThresholdFilterType::New();
    BinaryFilter->SetInput(result);

    //const short  outsidevalue = 0;
    //const short  insidevalue = 1000;

    BinaryFilter->SetOutsideValue(outsidevalue);
    BinaryFilter->SetInsideValue(insidevalue);

    //const short lowerThreshold = -800;
    //const short upperThreshold = 10000;

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
    //std::cout << labelMap->GetNumberOfLabelObjects() << " labels." << std::endl;
    //return ImageType::New();

    LabelType maxVolLabel = 0;
    itk::SizeValueType maxVol = 0;

    for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n)
    {
        ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);
        //std::cout << "maxVol = " << labelObject->GetNumberOfPixels() << std::endl;
        if (labelObject->GetNumberOfPixels() > maxVol)
        {
            maxVolLabel = labelObject->GetLabel();
            maxVol = labelObject->GetNumberOfPixels();
        }
    }
    //std::cout << "maxVolLabel = " << maxVolLabel << "," << "maxVol = " << maxVol << std::endl;

    ImageType::RegionType inputRegion = result->GetLargestPossibleRegion();
    ImageType::SizeType dims = inputRegion.GetSize();
    LabelType* labelImageImBuffer = connected->GetOutput()->GetBufferPointer();

    ImageType::PixelType* imBuffer = result->GetBufferPointer();
    //vtkSmartPointer<vtkPoints> points = getSplinePoints();
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

vtkSmartPointer<vtkPoints> AutomaticSegmentation::getSplinePoints()
{
    double p0[3] = { 0.0, -150.0, 0.0 };
    double p1[3] = { 100.0, -150.0, 0.0 };
    double p2[3] = { 100.0, 50.0, 0.0 };
    double p3[3] = { 0.0, 50.0, 0.0 };
    double p4[3] = { 0.0, -150.0, 0.0 };

    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(p0);
    points->InsertNextPoint(p1);
    points->InsertNextPoint(p2);
    points->InsertNextPoint(p3);
    points->InsertNextPoint(p4);

    vtkSmartPointer<vtkParametricSpline> spline =
        vtkSmartPointer<vtkParametricSpline>::New();
    spline->SetPoints(points);

    // Method1
#if(0)
    vtkSmartPointer<vtkPoints> outPoints =
        vtkSmartPointer<vtkPoints>::New();
    int pointsCnt = 11;
    double step = 1.0 / (pointsCnt - 1);
    for (double i = 0; i <= 1; i = i + step)
    {
        double u[] = { i, 0, 0 };
        double p[3];
        spline->Evaluate(u, p, NULL);
        outPoints->InsertNextPoint(p);
    }
#else
    // Method2
#define DIS(p1, p2) (sqrt((p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]))) 
    vtkSmartPointer<vtkPoints> outPoints =
        vtkSmartPointer<vtkPoints>::New();

    double totalLen = DIS(p0, p1) + DIS(p1, p2) + DIS(p2, p3) + DIS(p3, p4);
    double start = 0;
    double step = .2;

    while (start <= totalLen)
    {
        double u[] = { start / totalLen, 0, start / totalLen };
        double p[3];
        spline->Evaluate(u, p, NULL);
        outPoints->InsertNextPoint(p);
        start += step;
    }
#endif
    return outPoints;
#if (0)
    vtkSmartPointer<vtkParametricSpline> outSpline =
        vtkSmartPointer<vtkParametricSpline>::New();
    outSpline->SetPoints(outPoints);

    vtkSmartPointer<vtkParametricFunctionSource> functionSource =
        vtkSmartPointer<vtkParametricFunctionSource>::New();
    functionSource->SetParametricFunction(outSpline);
    functionSource->Update();

    vtkSmartPointer<vtkPolyDataMapper> splineMapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    splineMapper->SetInputConnection(functionSource->GetOutputPort());

    vtkSmartPointer<vtkActor> splineActor =
        vtkSmartPointer<vtkActor>::New();
    splineActor->SetMapper(splineMapper);

    vtkSmartPointer<vtkSphereSource> sphereSource =
        vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetPhiResolution(21);
    sphereSource->SetThetaResolution(21);
    sphereSource->SetRadius(.05);

    vtkSmartPointer<vtkPolyData> splinePointsData =
        vtkSmartPointer<vtkPolyData>::New();
    splinePointsData->SetPoints(outPoints);

    vtkSmartPointer<vtkGlyph3DMapper> splinePointsMapper =
        vtkSmartPointer<vtkGlyph3DMapper>::New();
    splinePointsMapper->SetInputData(splinePointsData);
    splinePointsMapper->SetSourceConnection(sphereSource->GetOutputPort());

    vtkSmartPointer<vtkActor> pointsActor =
        vtkSmartPointer<vtkActor>::New();
    pointsActor->SetMapper(splinePointsMapper);
    pointsActor->GetProperty()->SetColor(1, 0, 0);

    vtkSmartPointer<vtkRenderer> renderer =
        vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow =
        vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetSize(600, 600);
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
        vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
    ///////////////////////////////////////////

    renderer->AddActor(splineActor);
    renderer->AddActor(pointsActor);

    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
#endif
}


void
AutomaticSegmentation::DeepCopy(ImageType::Pointer input, ImageType::Pointer output)
{
    output->SetRegions(input->GetLargestPossibleRegion());
    output->Allocate();

    itk::ImageRegionConstIterator<ImageType> inputIterator(input, input->GetLargestPossibleRegion());
    itk::ImageRegionIterator<ImageType>      outputIterator(output, output->GetLargestPossibleRegion());

    while (!inputIterator.IsAtEnd())
    {
        outputIterator.Set(inputIterator.Get());
        ++inputIterator;
        ++outputIterator;
    }
}

SegmentImageType::Pointer AutomaticSegmentation::BinaryErodeImage3d(SegmentImageType* image, short r) const
{
    typedef itk::BinaryBallStructuringElement< SegmentImageType::PixelType, 3  > StructuringElementType;
    StructuringElementType  structuringElement;
    structuringElement.SetRadius(r);
    structuringElement.CreateStructuringElement();

    typedef itk::BinaryErodeImageFilter <
        SegmentImageType,
        SegmentImageType,
        StructuringElementType >  ErodeFilterType;

    ErodeFilterType::Pointer  binaryErode = ErodeFilterType::New();
    binaryErode->SetInput(image);
    binaryErode->SetKernel(structuringElement);
    binaryErode->SetErodeValue(1);
    binaryErode->Update();
    return binaryErode->GetOutput();
}


SegmentImageType::Pointer  AutomaticSegmentation::BinaryDilateImage3d(SegmentImageType* image, short r) const
{
    typedef itk::BinaryBallStructuringElement< SegmentImageType::PixelType, 3  > StructuringElementType;
    StructuringElementType  structuringElement;
    structuringElement.SetRadius(r);
    structuringElement.CreateStructuringElement();

    typedef itk::BinaryDilateImageFilter <
        SegmentImageType,
        SegmentImageType,
        StructuringElementType >  DilateFilterType;

    DilateFilterType::Pointer binaryDilate = DilateFilterType::New();
    binaryDilate->SetInput(image);
    binaryDilate->SetKernel(structuringElement);
    binaryDilate->SetDilateValue(1);
    binaryDilate->Update();

    //vis3d<ImageType>(binaryDilate->GetOutput());
    return binaryDilate->GetOutput();
}

SegmentImageType::Pointer AutomaticSegmentation::castImage(const ImageType::Pointer input) const
{
    constexpr unsigned int Dimension = 3;

    using InputImageType = itk::Image<PixelType, Dimension>;
    using OutputImageType = itk::Image<SegmentPixelType, Dimension>;

    using FilterType = itk::CastImageFilter<InputImageType, OutputImageType>;
    FilterType::Pointer filter = FilterType::New();

    filter->SetInput(input);
    filter->Update();
    return filter->GetOutput();
}
/*
SegmentImageType::Pointer AutomaticSegmentation::copyIm_byMark(ImageType* image, ImageType* binary_Im, short forePix)
{
    ImageType::RegionType inputRegion = image->GetLargestPossibleRegion();
    ImageType::SizeType dims = inputRegion.GetSize();
    ImageType::PixelType* imBuffer = image->GetBufferPointer();

    ImageType::PixelType* binaryImBuffer = binary_Im->GetBufferPointer();

    ImageType::PixelType a, b, c;
    bool act = false;
    a = forePix / 2;
    c = a;

    for (int k = 0; k < dims[2]; k++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int i = 0; i < dims[0]; i++)
            {
                int index = k * dims[1] * dims[0] + j * dims[0] + i;

                if (binaryImBuffer[index] == forePix)
                {
                    binaryImBuffer[index] = 1;
                }
                else
                {
                    binaryImBuffer[index] = 0;
                }

            }
        }
    }

    SegmentImageType::Pointer newImg = castImage(binary_Im);

    return newImg;
}
*/

ImageType::Pointer AutomaticSegmentation::OtsuMultipleThresholds(const ImageType* image, short &threshold)
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
    //std::cout << "filter->GetThresholds(): ";
    //std::cout << thresholds.size() << std::endl;
    //short  threshold = 300;
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

    //std::cout << threshold << std::endl;
    //threshold = 62;
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

/*
void AutomaticSegmentation::getBoneParts(SegmentImageType* image, SegmentImageType* tibiaImage, SegmentImageType* kneeCap)
{
    const unsigned int Dimension = 3;
    typedef uint8_t                                LabelType;
    typedef itk::Image< LabelType, Dimension >            OutputImageType;
    typedef itk::ShapeLabelObject< LabelType, Dimension > ShapeLabelObjectType;
    typedef itk::LabelMap< ShapeLabelObjectType >         LabelMapType;

    typedef itk::ConnectedComponentImageFilter <SegmentImageType, OutputImageType > ConnectedComponentImageFilterType;

    ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
    connected->SetInput(image);
    connected->Update();

    typedef itk::LabelImageToShapeLabelMapFilter< OutputImageType, LabelMapType> I2LType;
    I2LType::Pointer i2l = I2LType::New();
    i2l->SetInput(connected->GetOutput());
    i2l->SetComputePerimeter(true);
    i2l->Update();
    LabelMapType *labelMap = i2l->GetOutput();

    LabelType femurLabel = 0;
    LabelType tibiaLabel = 0;
    LabelType kneeLabel = 0;

    std::vector<LabelsInfo> LabelInfoVector;

    for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n)
    {
        ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);
        LabelInfoVector.push_back(LabelsInfo(labelObject->GetLabel(), labelObject->GetNumberOfPixels(), labelObject->GetRoundness(), labelObject->GetElongation(), labelObject->GetRegion()));
    }

    std::sort(LabelInfoVector.begin(), LabelInfoVector.end(), SortByPixels());

    ImageType::RegionType kneeRegion;

    if (LabelInfoVector.size() >= 3)
    {
        femurLabel = LabelInfoVector[0].MaxVolLabel;
        tibiaLabel = LabelInfoVector[1].MaxVolLabel;
        kneeLabel = LabelInfoVector[2].MaxVolLabel;
    }
    else
    {
        throw SegmentationException("It has not been possible to extract the bones of the femur, tibia and patella.");
    }

    SegmentImageType::RegionType inputRegion = image->GetLargestPossibleRegion();
    SegmentImageType::SizeType dims = inputRegion.GetSize();
    LabelType* labelImageImBuffer = connected->GetOutput()->GetBufferPointer();
    SegmentImageType::PixelType* imBuffer = image->GetBufferPointer();
    SegmentImageType::PixelType* tibiaImageBuffer = tibiaImage->GetBufferPointer();
    SegmentImageType::PixelType* kneeImageBuffer = kneeCap->GetBufferPointer();

    for (int k = 0; k < dims[2]; k++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int i = 0; i < dims[0]; i++)
            {
                int index = k * dims[1] * dims[0] + j * dims[0] + i;

                if (labelImageImBuffer[index] != femurLabel)
                {
                    imBuffer[index] = 0;
                }

                if (labelImageImBuffer[index] != tibiaLabel)
                {
                    tibiaImageBuffer[index] = 0;
                }

                if (labelImageImBuffer[index] != kneeLabel)
                {
                    kneeImageBuffer[index] = 0;
                }

            }
        }
    }
}

*/

SegmentImageType::Pointer AutomaticSegmentation::GetDigitallyReconstructedRadiograph(float pRotationX, float pRotationY, float pRotationZ, int pOutputSize, bool pVerbose) const
{
    // CT volume rotation around isocenter along x,y,z axis in degrees
    float rx = pRotationX;
    float ry = pRotationY;
    float rz = pRotationZ;

    // Translation parameter of the isocenter in mm
    float tx = 0.;
    float ty = 0.;
    float tz = 0.;

    // The pixel indices of the isocenter
    float cx = 0.;
    float cy = 0.;
    float cz = 0.;

    // Source to isocenter distance in mm
    float sid = 1000.;

    // Default pixel spacing in the iso-center plane in mm
    float sx = 1.;
    float sy = 1.;

    // Size of the output image in number of pixels
    int dx = pOutputSize;
    int dy = pOutputSize;

    // The central axis positions of the 2D images in continuous indices
    float o2Dx = 0;
    float o2Dy = 0;

    double threshold = 0;
    
    // Creation of a \code{ResampleImageFilter} enables coordinates for
    // each of the pixels in the DRR image to be generated. These
    // coordinates are used by the \code{RayCastInterpolateImageFunction}
    // to determine the equation of each corresponding ray which is cast
    // through the input volume.

    typedef itk::ResampleImageFilter<ImageType, ImageType > FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(leg_image_);
    filter->SetDefaultPixelValue(0);

    // An Euler transformation is defined to position the input volume.
    // The \code{ResampleImageFilter} uses this transform to position the
    // output DRR image for the desired view.

    typedef itk::CenteredEuler3DTransform< double >  TransformType;
    TransformType::Pointer transform = TransformType::New();
    transform->SetComputeZYX(true);
    TransformType::OutputVectorType translation;
    translation[0] = tx;
    translation[1] = ty;
    translation[2] = tz;
    // constant for converting degrees into radians
    const double dtr = (std::atan(1.0) * 4.0) / 180.0;
    transform->SetTranslation(translation);
    transform->SetRotation(dtr*rx, dtr*ry, dtr*rz);

    ImageType::PointType   imOrigin = leg_image_->GetOrigin();
    ImageType::SpacingType imRes = leg_image_->GetSpacing();

    typedef ImageType::RegionType     InputImageRegionType;
    typedef InputImageRegionType::SizeType InputImageSizeType;
    InputImageRegionType imRegion = leg_image_->GetBufferedRegion();
    InputImageSizeType   imSize = imRegion.GetSize();

    imOrigin[0] += imRes[0] * static_cast<double>(imSize[0]) / 2.0;
    imOrigin[1] += imRes[1] * static_cast<double>(imSize[1]) / 2.0;
    imOrigin[2] += imRes[2] * static_cast<double>(imSize[2]) / 2.0;

    TransformType::InputPointType center;
    center[0] = cx + imOrigin[0];
    center[1] = cy + imOrigin[1];
    center[2] = cz + imOrigin[2];

    transform->SetCenter(center);
    if (pVerbose)
    {
        std::cout << "Image size: "
            << imSize[0] << ", " << imSize[1] << ", " << imSize[2]
            << std::endl << "   resolution: "
            << imRes[0] << ", " << imRes[1] << ", " << imRes[2]
            << std::endl << "   origin: "
            << imOrigin[0] << ", " << imOrigin[1] << ", " << imOrigin[2]
            << std::endl << "   center: "
            << center[0] << ", " << center[1] << ", " << center[2]
            << std::endl << "Transform: " << transform << std::endl;
    }

    // The \code{RayCastInterpolateImageFunction} is instantiated and passed the transform
    // object. The \code{RayCastInterpolateImageFunction} uses this
    // transform to reposition the x-ray source such that the DRR image
    // and x-ray source move as one around the input volume. This coupling
    // mimics the rigid geometry of the x-ray gantry.

    typedef itk::RayCastInterpolateImageFunction<ImageType, double>
        InterpolatorType;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    interpolator->SetTransform(transform);

    // We can then specify a threshold above which the volume's
    // intensities will be integrated.

    interpolator->SetThreshold(threshold);

    // The ray-cast interpolator needs to know the initial position of the
    // ray source or focal point. In this example we place the input
    // volume at the origin and halfway between the ray source and the
    // screen. The distance between the ray source and the screen
    // is the "source to image distance" \code{sid} and is specified by
    // the user.

    InterpolatorType::InputPointType focalpoint;
    focalpoint[0] = imOrigin[0];
    focalpoint[1] = imOrigin[1];
    focalpoint[2] = imOrigin[2] - sid / 2.;
    interpolator->SetFocalPoint(focalpoint);
   
    // Having initialised the interpolator we pass the object to the
    // resample filter.

    if (pVerbose)
    {
        interpolator->Print(std::cout);
    }

    filter->SetInterpolator(interpolator);
    filter->SetTransform(transform);

    // The size and resolution of the output DRR image is specified via the
    // resample filter.

    ImageType::SizeType   size;
    size[0] = dx;  // number of pixels along X of the 2D DRR image
    size[1] = dy;  // number of pixels along Y of the 2D DRR image
    size[2] = 1;   // only one slice
    filter->SetSize(size);
    ImageType::SpacingType spacing;
    spacing[0] = sx;  // pixel spacing along X of the 2D DRR image [mm]
    spacing[1] = sy;  // pixel spacing along Y of the 2D DRR image [mm]
    spacing[2] = 1.0; // slice thickness of the 2D DRR image [mm]
    filter->SetOutputSpacing(spacing);
    // Software Guide : EndCodeSnippet
    if (pVerbose)
    {
        std::cout << "Output image size: "
            << size[0] << ", "
            << size[1] << ", "
            << size[2] << std::endl;
        std::cout << "Output image spacing: "
            << spacing[0] << ", "
            << spacing[1] << ", "
            << spacing[2] << std::endl;
    }

    // In addition the position of the DRR is specified. The default
    // position of the input volume, prior to its transformation is
    // half-way between the ray source and screen and unless specified
    // otherwise the normal from the "screen" to the ray source passes
    // directly through the centre of the DRR.

    double origin[3];
    origin[0] = imOrigin[0] + o2Dx - sx * ((double)dx - 1.) / 2.;
    origin[1] = imOrigin[1] + o2Dy - sy * ((double)dy - 1.) / 2.;
    origin[2] = imOrigin[2] + sid / 2.;
    filter->SetOutputOrigin(origin);
    // Software Guide : EndCodeSnippet
    if (pVerbose)
    {
        std::cout << "Output image origin: "
            << origin[0] << ", "
            << origin[1] << ", "
            << origin[2] << std::endl;
    }

    filter->Update();
    typedef itk::RescaleIntensityImageFilter<
        ImageType, SegmentImageType > RescaleFilterType;
    RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    rescaler->SetOutputMinimum(0);
    rescaler->SetOutputMaximum(255);
    rescaler->SetInput(filter->GetOutput());
    rescaler->Update();
    return rescaler->GetOutput();
}

/*

ImageType::Pointer AutomaticSegmentation::ReadSeries(const std::string& dirName, const std::string& seriesName)
{
    using PixelType = signed short;

    using NamesGeneratorType = itk::GDCMSeriesFileNames;
    NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

    nameGenerator->SetUseSeriesDetails(true);
    nameGenerator->AddSeriesRestriction("0008|0021");
    nameGenerator->SetGlobalWarningDisplay(false);
    nameGenerator->SetDirectory(dirName);

    ImageType::Pointer res = nullptr;

    try {
        using SeriesIdContainer = std::vector<std::string>;
        const SeriesIdContainer& seriesUID = nameGenerator->GetSeriesUIDs();
        auto seriesItr = seriesUID.begin();
        auto seriesEnd = seriesUID.end();

        if (seriesItr != seriesEnd) {
            //std::cout << "The directory: ";
            //std::cout << dirName << std::endl;
            //std::cout << "Contains the following DICOM Series: ";
            //std::cout << std::endl;
        }
        else {
            //std::cout << "No DICOMs in: " << dirName << std::endl;
            return nullptr;
        }
        while (seriesItr != seriesEnd) {
            //std::cout << seriesItr->c_str() << std::endl;
            ++seriesItr;
        }
        seriesItr = seriesUID.begin();
        while (seriesItr != seriesUID.end()) {
            //std::cout << "---------------------------------------------" << std::endl;
            std::string seriesIdentifier;
            if (!seriesName.empty()) // If seriesIdentifier given convert only that
            {
                seriesIdentifier = seriesName;
                seriesItr = seriesUID.end();
            }
            else // otherwise convert everything
            {
                seriesIdentifier = seriesItr->c_str();
                seriesItr++;
            }
            //std::cout << "\nReading: ";
            //std::cout << seriesIdentifier << std::endl;
            using FileNamesContainer = std::vector<std::string>;
            FileNamesContainer fileNames = nameGenerator->GetFileNames(seriesIdentifier);
            if (fileNames.empty()) {
                continue;
            }

            using ReaderType = itk::ImageSeriesReader<ImageType>;
            ReaderType::Pointer reader = ReaderType::New();
            using ImageIOType = itk::GDCMImageIO;
            ImageIOType::Pointer dicomIO = ImageIOType::New();
            reader->SetImageIO(dicomIO);
            reader->SetFileNames(fileNames);
            // reader->ForceOrthogonalDirectionOff(); // properly read CTs with gantry tilt

            reader->Update();
            // Read only first series
            //typename ImageType::Pointer res = reader->GetOutput();
            res = reader->GetOutput();
        }
    }
    catch (itk::ExceptionObject& ex) {
        std::cout << ex << std::endl;
        return res;

    }

    return res;
}



void AutomaticSegmentation::vis3d(ImageType* image)
{
    auto imageToVtkImage = itk::ImageToVTKImageFilter<ImageType>::New();
    imageToVtkImage->SetInput(image);
    imageToVtkImage->Update();
    //////////////////////////////////////
    vtkNew<vtkImageToPolyDataFilter> i2pd;
    i2pd->SetInputData(imageToVtkImage->GetOutput());

    vtkSmartPointer<vtkPolyDataMapper> outlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    outlineMapper->SetInputData(i2pd->GetOutput());

    vtkSmartPointer<vtkActor> outlineActor = vtkSmartPointer<vtkActor>::New();
    outlineActor->SetMapper(outlineMapper);

    vtkSmartPointer<vtkRenderer> ren[4];
    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
    renWin->SetMultiSamples(0);

    for (int i = 0; i < 4; i++)
    {
        ren[i] = vtkSmartPointer<vtkRenderer>::New();
        renWin->AddRenderer(ren[i]);
    }

    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);

    vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(0.005);

    vtkSmartPointer<vtkProperty> ipwProp = vtkSmartPointer<vtkProperty>::New();

    vtkSmartPointer<vtkImagePlaneWidget> planeWidget[3];
    int imageDims[3];
    imageToVtkImage->GetOutput()->GetDimensions(imageDims);

    for (int i = 0; i < 3; i++)
    {
        planeWidget[i] = vtkSmartPointer<vtkImagePlaneWidget>::New();
        planeWidget[i]->SetInteractor(iren);
        planeWidget[i]->SetPicker(picker);
        planeWidget[i]->RestrictPlaneToVolumeOn();
        double color[3] = { 0, 0, 0 };
        color[i] = 1;
        planeWidget[i]->GetPlaneProperty()->SetColor(color);
        planeWidget[i]->SetTexturePlaneProperty(ipwProp);
        planeWidget[i]->TextureInterpolateOff();
        planeWidget[i]->SetResliceInterpolateToLinear();
        planeWidget[i]->SetInputData(imageToVtkImage->GetOutput());
        planeWidget[i]->SetPlaneOrientation(i);
        planeWidget[i]->SetSliceIndex(imageDims[i] / 2);
        planeWidget[i]->DisplayTextOn();
        planeWidget[i]->SetDefaultRenderer(ren[3]);
        planeWidget[i]->SetWindowLevel(1358, -27);
        planeWidget[i]->On();
        planeWidget[i]->InteractionOn();
    }

    planeWidget[1]->SetLookupTable(planeWidget[0]->GetLookupTable());
    planeWidget[2]->SetLookupTable(planeWidget[0]->GetLookupTable());

    vtkSmartPointer<vtkResliceCursorCallback> cbk = vtkSmartPointer<vtkResliceCursorCallback>::New();
    vtkSmartPointer<vtkResliceCursor> resliceCursor = vtkSmartPointer<vtkResliceCursor >::New();
    resliceCursor->SetCenter(imageToVtkImage->GetOutput()->GetCenter());
    resliceCursor->SetThickMode(0);
    resliceCursor->SetThickness(10, 10, 10);
    resliceCursor->SetImage(imageToVtkImage->GetOutput());

    vtkSmartPointer< vtkResliceCursorWidget > resliceCursorWidget[3];
    vtkSmartPointer< vtkResliceCursorLineRepresentation > resliceCursorRep[3];

    double viewUp[3][3] = { { 0, 0, -1 }, { 0, 0, 1 }, { 0, 1, 0 } };
    for (int i = 0; i < 3; i++)
    {
        resliceCursorWidget[i] = vtkSmartPointer< vtkResliceCursorWidget >::New();
        resliceCursorWidget[i]->SetInteractor(iren);

        resliceCursorRep[i] = vtkSmartPointer< vtkResliceCursorLineRepresentation >::New();
        resliceCursorWidget[i]->SetRepresentation(resliceCursorRep[i]);
        resliceCursorRep[i]->GetResliceCursorActor()->GetCursorAlgorithm()->SetResliceCursor(resliceCursor);
        resliceCursorRep[i]->GetResliceCursorActor()->GetCursorAlgorithm()->SetReslicePlaneNormal(i);

        const double minVal = imageToVtkImage->GetOutput()->GetScalarRange()[0];
        if (vtkImageReslice *reslice = vtkImageReslice::SafeDownCast(resliceCursorRep[i]->GetReslice()))
        {
            reslice->SetBackgroundColor(minVal, minVal, minVal, minVal);
        }

        resliceCursorWidget[i]->SetDefaultRenderer(ren[i]);
        resliceCursorWidget[i]->SetEnabled(1);

        ren[i]->GetActiveCamera()->SetFocalPoint(0, 0, 0);
        double camPos[3] = { 0, 0, 0 };
        camPos[i] = 1;
        ren[i]->GetActiveCamera()->SetPosition(camPos);
        ren[i]->GetActiveCamera()->ParallelProjectionOn();
        ren[i]->GetActiveCamera()->SetViewUp(viewUp[i][0], viewUp[i][1], viewUp[i][2]);
        ren[i]->ResetCamera();
        cbk->IPW[i] = planeWidget[i];
        cbk->RCW[i] = resliceCursorWidget[i];
        resliceCursorWidget[i]->AddObserver(vtkResliceCursorWidget::ResliceAxesChangedEvent, cbk);
        double range[2];
        imageToVtkImage->GetOutput()->GetScalarRange(range);
        resliceCursorRep[i]->SetWindowLevel(range[1] - range[0], (range[0] + range[1]) / 2.0);
        planeWidget[i]->SetWindowLevel(range[1] - range[0], (range[0] + range[1]) / 2.0);
        resliceCursorRep[i]->SetLookupTable(resliceCursorRep[0]->GetLookupTable());
        planeWidget[i]->GetColorMap()->SetLookupTable(resliceCursorRep[0]->GetLookupTable());
    }

    ren[0]->SetBackground(0.3, 0.1, 0.1);
    ren[1]->SetBackground(0.1, 0.3, 0.1);
    ren[2]->SetBackground(0.1, 0.1, 0.3);
    ren[3]->AddActor(outlineActor);
    ren[3]->SetBackground(0.1, 0.1, 0.1);
    renWin->SetSize(600, 600);

    ren[0]->SetViewport(0, 0, 0.5, 0.5);
    ren[1]->SetViewport(0.5, 0, 1, 0.5);
    ren[2]->SetViewport(0, 0.5, 0.5, 1);
    ren[3]->SetViewport(0.5, 0.5, 1, 1);
    renWin->Render();

    ren[3]->GetActiveCamera()->Elevation(110);
    ren[3]->GetActiveCamera()->SetViewUp(0, 0, -1);
    ren[3]->GetActiveCamera()->Azimuth(45);
    ren[3]->GetActiveCamera()->Dolly(1.15);
    ren[3]->ResetCameraClippingRange();

    vtkSmartPointer< vtkInteractorStyleImage > style = vtkSmartPointer< vtkInteractorStyleImage >::New();
    iren->SetInteractorStyle(style);
    iren->Initialize();
    iren->Start();
    //return EXIT_SUCCESS;
}*/