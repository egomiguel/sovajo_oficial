#ifndef MANUAL_SEGMENTATION_H
#define MANUAL_SEGMENTATION_H

#include "Types.hpp"
#include "tka_segmentation_export.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"

namespace TKA
{
	namespace SEGMENTATION
	{

		const int8_t BinaryBackgroundPixel = 0;
		const int8_t BinaryNoBackgroundPixel = 1;

		class TKA_SEGMENTATION_EXPORT ManualSegmentation
		{
		public:

			enum ViewPlane {
				kSagital_X,
				kCoronal_Y,
				kAxial_Z
			};

			ManualSegmentation(const ImageType::Pointer pImage);

			ImageType::Pointer ApplyThresholds(short pMinthreshold, short pMaxthreshold, bool pDenoiseBefore = true);

			static vtkSmartPointer<vtkPolyData> BinaryImageToPolyData(const SegmentImageType::Pointer pBinaryImage);
			/*
				pAnatomicalPlane: Cutting plane in its anatomical version, that is sagittal, coronal or axial.

				pObliqueVector: Vector of the cutting plane in its original version, that is, inclined.
			*/
			static ImageType::Pointer MakeRotation(const ImageType::Pointer& pImageIn, const ViewPlane& pAnatomicalPlane, const itk::Vector<double, 3>& pObliqueVector, itk::Matrix< double, 3, 3 >& pRotationOut, ImageType::PixelType pDefaultPixelValue = ImageTypeMin);

			static SegmentImageType::Pointer RestoreImageSegmentation(const ImageType::Pointer& pSegmentationIn, const itk::Matrix< double, 3, 3 >& pRotationIn);

			/*
				Make a flood fill in closed contours. The outer area of the contours (auter background) must be greater than each inner area.
			*/
			static ImageType::Pointer FloodFillBinaryImageSlice(const ImageType::Pointer& pImage, int8_t pNoBackgroundPixel = BinaryNoBackgroundPixel);

			static ImageType2D::Pointer FloodFillBinaryImageSlice(const ImageType2D::Pointer& pImage, int8_t pNoBackgroundPixel = BinaryNoBackgroundPixel);

		private:

			ImageType::Pointer mLegImage;

			void DeepCopy(const ImageType::Pointer pInput, ImageType::Pointer pOutput) const;

			static SegmentImageType::Pointer CastImage(const ImageType::Pointer pInput);

			ImageType::Pointer GradientAnisotropicDiffusion(const ImageType::Pointer pImage) const;

			template<typename ImageType>
			static typename ImageType::Pointer CloneImage(const typename ImageType::Pointer pInput);

			template<typename ImageType>
			static typename ImageType::Pointer FillHole(typename ImageType::Pointer pImage, typename ImageType::PixelType pValue);

			template<typename ImageType>
			static typename ImageType::Pointer RotateImage(typename ImageType::Pointer pImage, const itk::Matrix< double, 3, 3 >& pRotationIn, typename ImageType::PixelType pDefaultPixelValue);

			void Get2DSlice(const ImageType::Pointer& pImageIn, ImageType::Pointer& pImageOut, int pSliceNumber, const ViewPlane& pPlane) const;
		};
	}
}

#endif