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
#include "itkImageToVTKImageFilter.hpp"
#include "itkRelabelComponentImageFilter.h"

#include "SegmentationException.hpp"
#include "BallsManualSegmentation.hpp"

using namespace SPINE::SEGMENTATION;


BallsManualSegmentation::BallsManualSegmentation(const ImageType::Pointer pImage)
{
	mImage = pImage;
}

SegmentImageType::Pointer BallsManualSegmentation::getSegmentBallsAndCentroids(short pMinthreshold, short pMaxthreshold, short pMinSize, std::vector<itk::Point<double, 3>>& pCentroids) const
{
	auto thresholdFilter = itk::BinaryThresholdImageFilter<ImageType, ImageType>::New();
	thresholdFilter->SetInput(mImage);
	thresholdFilter->SetUpperThreshold(pMaxthreshold);
	thresholdFilter->SetLowerThreshold(pMinthreshold);
	thresholdFilter->SetInsideValue(BinaryNoBackgroundPixel);
	thresholdFilter->SetOutsideValue(BinaryBackgroundPixel);

	using LabelImageType = itk::Image<PixelType, 3>;
	auto connectedFilter = itk::ConnectedComponentImageFilter<ImageType, LabelImageType>::New();
	connectedFilter->SetInput(thresholdFilter->GetOutput());

	auto relabelFilter = itk::RelabelComponentImageFilter<LabelImageType, LabelImageType>::New();
	relabelFilter->SetInput(connectedFilter->GetOutput());
	relabelFilter->SetMinimumObjectSize(pMinSize);

	using LabelMapType = itk::LabelMap<itk::ShapeLabelObject<PixelType, 3>>;
	auto labelMapFilter = itk::LabelImageToShapeLabelMapFilter<LabelImageType, LabelMapType>::New();
	labelMapFilter->SetInput(relabelFilter->GetOutput());
	labelMapFilter->Update();

	LabelMapType *labelMap = labelMapFilter->GetOutput();
	std::vector<unsigned long> sizes;
	for (unsigned int i = 0; i < labelMap->GetNumberOfLabelObjects(); ++i) {
		auto labelObject = labelMap->GetNthLabelObject(i);
		sizes.push_back(labelObject->GetNumberOfPixels());
	}

	std::sort(sizes.begin(), sizes.end());
	unsigned long medianSize = sizes[sizes.size() / 2];
	//unsigned long rangeMin = static_cast<unsigned long>(medianSize * 0.8);
	unsigned long rangeMax = static_cast<unsigned long>(medianSize * 2);

	auto filteredLabels = LabelImageType::New();
	filteredLabels->SetRegions(relabelFilter->GetOutput()->GetLargestPossibleRegion());
	filteredLabels->Allocate();
	filteredLabels->FillBuffer(0);

	for (unsigned int i = 0; i < labelMap->GetNumberOfLabelObjects(); ++i) {
		auto labelObject = labelMap->GetNthLabelObject(i);
		unsigned long size = labelObject->GetNumberOfPixels();

		if (size <= rangeMax) {
			auto boundingBox = labelObject->GetBoundingBox();
			pCentroids.push_back(labelObject->GetCentroid());

			itk::ImageRegionIterator<LabelImageType> it(filteredLabels, boundingBox);
			for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
				if (labelObject->HasIndex(it.GetIndex())) {
					it.Set(labelObject->GetLabel());
				}
			}
		}
	}

	auto maskFilter = itk::MaskImageFilter<ImageType, LabelImageType, ImageType>::New();
	maskFilter->SetInput(mImage);
	maskFilter->SetMaskImage(filteredLabels);
	return CastImage(maskFilter->GetOutput());

}

SegmentImageType::Pointer BallsManualSegmentation::CastImage(const ImageType::Pointer pInput)
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
