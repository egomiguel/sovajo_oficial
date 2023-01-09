#ifndef AUTOMATIC_SEGMENTATION_H
#define AUTOMATIC_SEGMENTATION_H

#include <string>
#include "ResliceCallbackVTK.hpp"
#include "Types.hpp"
#include "segmentation_export.h"

class SEGMENTATION_EXPORT AutomaticSegmentation
{
public:
    AutomaticSegmentation(ImageType::Pointer legImage);

    bool ExecuteSegmentation();

    SegmentImageType::Pointer GetFemoralSegment() const;

    SegmentImageType::Pointer GetTibiaSegment() const;

    SegmentImageType::Pointer GetKneeCapSegment() const;

    SegmentImageType::Pointer GetFibulaSegment() const;

    SegmentImageType::Pointer GetDigitallyReconstructedRadiograph(float pRotationX = 0, float pRotationY = 0, float pRotationZ = 0, int pOutputSize = 1024, bool pVerbose = false) const;

	//void vis3d(ImageType* image);

private:

    ImageType::Pointer leg_image_, leg_image_clean_;

    SegmentImageType::Pointer femur_segment_, tibia_segment_, knee_cap_, fibula_segment;

    void getSegmentationBones(const ImageType::Pointer image);

    SegmentImageType::Pointer closeBone(const SegmentImageType::Pointer image, short r) const;

	vtkSmartPointer<vtkPoints> getSplinePoints();

	ImageType::Pointer getOnlyLeg(const ImageType::Pointer image, PixelType outsidevalue, PixelType insidevalue,
		PixelType lowerThreshold, PixelType upperThreshold);

    double getSmoothParameter(const ImageType::Pointer image, PixelType outsidevalue, PixelType insidevalue,
        PixelType lowerThreshold, PixelType upperThreshold) const;

	ImageType::Pointer OtsuMultipleThresholds(const ImageType* image, short &threshold);

    //ImageType::Pointer getMarkers(const ImageType::Pointer image, ImageType::Pointer cleanImage) const;

	void DeepCopy(ImageType::Pointer input, ImageType::Pointer output);

    template<typename ImageType>
    typename ImageType::Pointer CloneImage(const typename ImageType::Pointer input) const;

	//void getBoneParts(SegmentImageType* image, SegmentImageType* tibiaImage, SegmentImageType* kneeCap);

    SegmentImageType::Pointer BinaryErodeImage3d(SegmentImageType* image, short r) const;

    SegmentImageType::Pointer BinaryDilateImage3d(SegmentImageType* image, short r) const;

    //SegmentImageType::Pointer copyIm_byMark(ImageType* image, ImageType* binary_Im, short forePix);

    SegmentImageType::Pointer castImage(const ImageType::Pointer input) const;

    SegmentImageType::Pointer WatershedSegmentation(ImageType::Pointer image, ImageType::Pointer marker) const;

    ImageType::Pointer GradientAnisotropicDiffusion(const ImageType::Pointer image);
};

#endif