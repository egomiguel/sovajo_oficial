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
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "SpineLine.hpp"
#include "itkGradientMagnitudeImageFilter.h"
#include <numeric>
#include <cmath>
#include "SegmentationException.hpp"
#include "SpineSegmentation.hpp"
#include <fstream>
#include <stdio.h>

using namespace SPINE::SEGMENTATION;

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


inline double getAngleBetweenVectors_Spine(const cv::Point3d& a, const cv::Point3d& b)
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


inline cv::Mat getRotateMatrix_Spine(const cv::Point3d& axis, double angle)
{
	cv::Mat rotationMatrix(3, 3, CV_64F);

	cv::Point3d normaliceAxis = axis / sqrt(axis.dot(axis));;

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

inline cv::Mat GetRotateLineToX_Spine(const cv::Point3d& vector)
{
	cv::Point3d normalXY(1, 0, 0);
	cv::Point3d rotationAxis = vector.cross(normalXY);
	rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
	double rotationAngle = getAngleBetweenVectors_Spine(normalXY, vector);
	cv::Mat rotation = getRotateMatrix_Spine(rotationAxis, rotationAngle);
	return rotation;
}

inline double getDistancePoints(const cv::Point3d& A, const cv::Point3d& B)
{
	auto diff = B - A;
	return sqrt(diff.dot(diff));
}


SpineSegmentation::SpineSegmentation(ImageType::Pointer pSpine)
{
    mSpine = pSpine;
}

cv::Point3d SpineSegmentation::getFixNormal(const cv::Point3d& V1, const cv::Point3d& V2, bool fixFirst) const
{
	double degreeThreshold = 19.;
	double angle_rad = getAngleBetweenVectors_Spine(V1, V2);
	double angleDegree = angle_rad * 180. / 3.14159;
	if (180. - angleDegree > degreeThreshold)
	{
		double diff = (180. - angleDegree) - degreeThreshold;
		diff = diff * 3.14159 / 180.;
		cv::Point3d axis, finalNormal;
		if (fixFirst == true)
		{
			axis = V2.cross(V1);
		}
		else
		{
			axis = V1.cross(V2);
		}

		auto rotation = getRotateMatrix_Spine(axis, diff);

		if (fixFirst == true)
		{
			cv::Mat pointMat(3, 1, CV_64F);
			pointMat.at <double>(0, 0) = V1.x;
			pointMat.at <double>(1, 0) = V1.y;
			pointMat.at <double>(2, 0) = V1.z;

			cv::Mat rotateV1 = rotation * pointMat;
			cv::Point3d V1_temp = cv::Point3d(rotateV1);
			V1_temp = V1_temp / sqrt(V1_temp.dot(V1_temp));
			return V1_temp;
		}
		else
		{
			cv::Mat pointMat(3, 1, CV_64F);
			pointMat.at <double>(0, 0) = V2.x;
			pointMat.at <double>(1, 0) = V2.y;
			pointMat.at <double>(2, 0) = V2.z;

			cv::Mat rotateV2 = rotation * pointMat;
			cv::Point3d V2_temp = cv::Point3d(rotateV2);
			V2_temp = V2_temp / sqrt(V2_temp.dot(V2_temp));
			return V2_temp;
		}
	}
	else
	{
		if (fixFirst == true)
		{
			return V1;
		}
		else
		{
			return V2;
		}
	}
}

/*
	The highest correlation is sought between the series of points in the center of the column. 
	This correlation defines the size of a sliding window within which the pixel with the lowest value is chosen.
*/

std::vector<SpineSegmentation::Plane> SpineSegmentation::getIntervertebralPlanes(const std::vector<ImageType::PointType>& centerPhysicalPoints, const ImageType::RegionType& region) const
{
	std::vector<SpineSegmentation::Plane> result;

	std::vector<double> pixels;

	std::vector<ImageType::PointType> morePoints;

	for (int i = 0; i < centerPhysicalPoints.size() - 1; i++)
	{
		auto currentPoint = cv::Point3d(centerPhysicalPoints[i][0], centerPhysicalPoints[i][1], centerPhysicalPoints[i][2]);
		auto nextPoint = cv::Point3d(centerPhysicalPoints[i + 1][0], centerPhysicalPoints[i + 1][1], centerPhysicalPoints[i + 1][2]);
		auto vector = nextPoint - currentPoint;
		vector = vector / sqrt(vector.dot(vector));

		ImageType::PointType physicalPoint1;
		physicalPoint1[0] = currentPoint.x;
		physicalPoint1[1] = currentPoint.y;
		physicalPoint1[2] = currentPoint.z;
		morePoints.push_back(physicalPoint1);

		auto tempPoint = currentPoint;
		while (getDistancePoints(tempPoint, nextPoint) > 1.3)
		{
			tempPoint = tempPoint + vector;
			ImageType::PointType physicalPoint2;
			physicalPoint2[0] = tempPoint.x;
			physicalPoint2[1] = tempPoint.y;
			physicalPoint2[2] = tempPoint.z;
			morePoints.push_back(physicalPoint2);
		}
	}

	morePoints.push_back(centerPhysicalPoints[centerPhysicalPoints.size() - 1]);

	for (auto& item : morePoints)
	{
		ImageType::IndexType indexPoint;
		if (mSpine->TransformPhysicalPointToIndex(item, indexPoint))
		{
			std::vector<double> values;
			values.push_back(getPixel(mSpine, indexPoint[0], indexPoint[1], indexPoint[2]));

			values.push_back(getPixel(mSpine, indexPoint[0] + 1, indexPoint[1], indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] + 1, indexPoint[1] + 1, indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] + 1, indexPoint[1] - 1, indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] + 1, indexPoint[1] + 2, indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] + 1, indexPoint[1] - 2, indexPoint[2]));

			values.push_back(getPixel(mSpine, indexPoint[0] + 2, indexPoint[1], indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] + 2, indexPoint[1] + 1, indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] + 2, indexPoint[1] - 1, indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] + 2, indexPoint[1] + 2, indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] + 2, indexPoint[1] - 2, indexPoint[2]));

			values.push_back(getPixel(mSpine, indexPoint[0] - 1, indexPoint[1], indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] - 1, indexPoint[1] + 1, indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] - 1, indexPoint[1] - 1, indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] - 1, indexPoint[1] + 2, indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] - 1, indexPoint[1] - 2, indexPoint[2]));

			values.push_back(getPixel(mSpine, indexPoint[0] - 2, indexPoint[1], indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] - 2, indexPoint[1] + 1, indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] - 2, indexPoint[1] - 1, indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] - 2, indexPoint[1] + 2, indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0] - 2, indexPoint[1] - 2, indexPoint[2]));

			values.push_back(getPixel(mSpine, indexPoint[0], indexPoint[1] + 1, indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0], indexPoint[1] - 1, indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0], indexPoint[1] + 2, indexPoint[2]));
			values.push_back(getPixel(mSpine, indexPoint[0], indexPoint[1] - 2, indexPoint[2]));

			values.push_back(getPixel(mSpine, indexPoint[0], indexPoint[1], indexPoint[2] + 1));
			values.push_back(getPixel(mSpine, indexPoint[0], indexPoint[1], indexPoint[2] - 1));
			values.push_back(getPixel(mSpine, indexPoint[0], indexPoint[1], indexPoint[2] + 2));
			values.push_back(getPixel(mSpine, indexPoint[0], indexPoint[1], indexPoint[2] - 2));

			double avg = calculateMean(values);
			pixels.push_back(avg);
		}
		else
		{
			throw SegmentationExceptionCode::PHYSICAL_POINT_TO_INDEX_OUTSIDE_OF_IMAGE;
		}
	}

	double bestCorr = -1;
	int bestCorrlag = -1;
	for (int lag = 20; lag <= 60; lag++)
	{
		if (lag >= pixels.size())
		{
			break;
		}
		double correlation = calculateCorrelation(pixels, lag);
		if (correlation > bestCorr)
		{
			bestCorr = correlation;
			bestCorrlag = lag;
		}
	}

	if (bestCorrlag == -1)
	{
		return result;
	}
	//std::cout << "Size: " << pixels.size() << ", Correlation: " << bestCorrlag << std::endl;

	int generalMin = getLowerValueAsPos(pixels, 0, pixels.size() - 1);
	std::set<int> allPos;
	allPos.insert(generalMin);

	int tempPos = generalMin;

	double lagOffset = bestCorrlag / 10;

	/*
		Searching for the lowest values ​​from the general minimum backwards and using the correlation 
		as the search size. The correlation is shifted by (correlation/2) and from there a window of 
		size equal to the correlation is searched.
	*/

	while (tempPos > 0)
	{
		int tEnd = tempPos - (bestCorrlag / 2);
		int tBegin = tEnd - bestCorrlag;

		//std::cout << "Inicio: " << tBegin << " tEnd: " << tEnd << "Total: "<< pixels.size() << std::endl;

		if (tBegin < 0 && (tempPos - bestCorrlag >= 0 || abs(tempPos - bestCorrlag) <= lagOffset))
		{
			tBegin = 0;
		}

		if (tBegin < 0)
		{
			break;
		}

		tempPos = getLowerValueAsPos(pixels, tBegin, tEnd);
		allPos.insert(tempPos);
	}

	/*
		Searching for the lowest values ​​from the general minimum forwards and using the correlation
		as the search size. The correlation is shifted by (correlation/2) and from there a window of
		size equal to the correlation is searched.
	*/

	tempPos = generalMin;

	while (tempPos < pixels.size())
	{
		int tBegin = tempPos + (bestCorrlag / 2);
		int tEnd = tBegin + bestCorrlag;

		//std::cout << "Inicio2: " << tBegin << " tEnd2: " << tEnd << "Total: " << pixels.size() << std::endl;
		
		if (tEnd >= pixels.size() && ((tempPos + bestCorrlag) - pixels.size()) <= lagOffset)
		{
			tEnd = pixels.size() - 1;
		}

		if (tEnd >= pixels.size())
		{
			break;
		}

		tempPos = getLowerValueAsPos(pixels, tBegin, tEnd);
		allPos.insert(tempPos);
	}

	/*
	auto it1 = allPos.begin();
	auto it2 = allPos.end();

	for ( ; it1 != it2; ++it1)
	{
		std::cout << *it1 << std::endl;
	}

	auto imageTemp = mSpine;

	ImageType::RegionType inputRegion = imageTemp->GetLargestPossibleRegion();
	ImageType::SizeType dims = inputRegion.GetSize();

	ImageType::PixelType* imBuffer = imageTemp->GetBufferPointer();

	it1 = allPos.begin();
	it2 = allPos.end();

	for (; it1 != it2; ++it1)
	{
		for (int k = 0; k < dims[2]; k++)
		{
			for (int j = 0; j < dims[1]; j++)
			{
				for (int i = 0; i < dims[0]; i++)
				{
					int index = k * dims[1] * dims[0] + j * dims[0] + i;

					auto ref = physicalPoints[*it1];

					ImageType::IndexType indexPoint;

					if (imageTemp->TransformPhysicalPointToIndex(ref, indexPoint))
					{
						cv::Point3d miPunto1(indexPoint[0], indexPoint[1], indexPoint[2]);
						cv::Point3d miPunto2(i, j, k);

						double dist = sqrt((miPunto2 - miPunto1).dot(miPunto2 - miPunto1));

						if (dist < 6)
						{
							imBuffer[index] = 1000;
						}
					}
				}
			}
		}
	}

	SaveImageTest<ImageType>(imageTemp, "lastImage");
	*/

	auto it1 = allPos.begin();
	auto it2 = allPos.end();

	std::vector<cv::Point3d> centers;

	for (; it1 != it2; ++it1)
	{
		auto temp = morePoints[*it1];
		centers.push_back(cv::Point3d(temp[0], temp[1], temp[2]));
		/*
		ImageType::IndexType indexPoint;
		if (mSpine->TransformPhysicalPointToIndex(temp, indexPoint))
		{
			std::cout << "Small: " << getPixel(mSpine, indexPoint[0], indexPoint[1], indexPoint[2]) << std::endl;
		}
		*/
	}

	cv::Point3d beforeNormal;

	for (int i = 0; i < centers.size(); i++)
	{
		cv::Point3d normalTemp, V1, V2;

		if (i == 0)
		{
			V1 = centers[i] - centers[i + 1];
			V1 = V1 / sqrt(V1.dot(V1));

			normalTemp = V1;
			//std::cout << normalTemp << std::endl;
			beforeNormal = normalTemp;
		}
		else
		{
			if (i < centers.size() - 1)
			{
				V2 = centers[i + 1] - centers[i];
				normalTemp = -getFixNormal(beforeNormal, V2, false);
				normalTemp = normalTemp / sqrt(normalTemp.dot(normalTemp));
				//std::cout << normalTemp << std::endl;
				beforeNormal = normalTemp;
			}
			else
			{
				normalTemp = beforeNormal;
			}
		}
		

		itk::Point<double, 3> normal;
		itk::Point<double, 3> center;

		normal[0] = normalTemp.x;
		normal[1] = normalTemp.y;
		normal[2] = normalTemp.z;

		center[0] = centers[i].x;
		center[1] = centers[i].y;
		center[2] = centers[i].z;

		SpineSegmentation::Plane obj;
		obj.normal = normal;
		obj.center = center;

		result.push_back(obj);
	}

	/*
	for (int i = 0; i < result.size() - 1; i++)
	{
		auto n1 = result[i].normal;
		auto n2 = result[i + 1].normal;

		cv::Point3d V1 = cv::Point3d(n1[0], n1[1], n1[2]);
		cv::Point3d V2 = -cv::Point3d(n2[0], n2[1], n2[2]);


		double angle = getAngleBetweenVectors_Spine(V1, V2);
		std::cout <<V1 << "; " <<V2 << " Angle: " << angle * 180 / 3.14 << std::endl;
		
	}
	*/
	
	return result;
}

double SpineSegmentation::getPixel(const ImageType::Pointer image, int x, int y, int z) const
{
	ImageType::IndexType indexPoint;
	indexPoint[0] = x;
	indexPoint[1] = y;
	indexPoint[2] = z;

	return image->GetPixel(indexPoint);
}

/*
ImageType::Pointer SpineSegmentation::OtsuMultipleThresholds(const ImageType::Pointer image, short &threshold) const
{
	auto smoothingFilter = itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType>::New();
	smoothingFilter->SetSigma(0.5);
	smoothingFilter->SetInput(image);
	smoothingFilter->Update();

	ImageType::Pointer smoothingImage = smoothingFilter->GetOutput(); 

	cv::Point3d A(-10.5, -7.9, -313.24);
	cv::Point3d B(-10.35, -9.66, -467.24);
	SpineLine centerLine = SpineLine::makeLineWithPoints(A, B);

	std::vector<cv::Point2d> allPoints;
	std::vector<double> pixels;
	double distance = sqrt((A - B).dot(A - B));
	std::cout << "Distance: " << distance << std::endl;
	double step = 1. / distance;
	for (float i = 0; i < 1 + step; i += step)
	{
		auto tempPoint = A + i * (B - A);
		ImageType::PointType physicalPoint;
		physicalPoint[0] = tempPoint.x;
		physicalPoint[1] = tempPoint.y;
		physicalPoint[2] = tempPoint.z;

		ImageType::IndexType indexPoint;
		if (image->TransformPhysicalPointToIndex(physicalPoint, indexPoint))
		{
			std::vector<float> values;
			values.push_back(getPixel(image, indexPoint[0], indexPoint[1], indexPoint[2]));

			values.push_back(getPixel(image, indexPoint[0] + 1, indexPoint[1], indexPoint[2]));
			values.push_back(getPixel(image, indexPoint[0] - 1, indexPoint[1], indexPoint[2]));

			values.push_back(getPixel(image, indexPoint[0], indexPoint[1] + 1, indexPoint[2]));
			values.push_back(getPixel(image, indexPoint[0], indexPoint[1] - 1, indexPoint[2]));

			values.push_back(getPixel(image, indexPoint[0], indexPoint[1], indexPoint[2] + 1));
			values.push_back(getPixel(image, indexPoint[0], indexPoint[1], indexPoint[2] - 1));

			int sum = 0; 
			for (int item : values) 
			{
				sum += item;
			}

			allPoints.push_back(cv::Point2d(i * distance, sum /float(values.size()) ));
			pixels.push_back(sum / float(values.size()));

		}
	}

	
	std::cout << "init" << std::endl;
	std::ofstream file("puntos.txt");

	// Itera sobre el vector y escribe cada punto en el fichero
	for (const auto& point : allPoints) {
		file << point.x << "," << point.y << "\n";
	}

	// Cierra el fichero
	file.close();

	for (int lag = 20; lag <= 60; lag++)
	{
		double correlation = calculateCorrelation(pixels, lag);
		std::cout << "Lag: " << lag << " , " << correlation << std::endl;
	}

	std::cout << "end" << std::endl;

    return smoothingImage;
}
*/


double SpineSegmentation::calculateMean(const std::vector<double>& values) const {
	return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

double SpineSegmentation::calculateCovariance(const std::vector<double>& x, const std::vector<double>& y, double meanX, double meanY) const {
	double covariance = 0.0;
	for (size_t i = 0; i < x.size(); ++i) {
		covariance += (x[i] - meanX) * (y[i] - meanY);
	}
	return covariance / x.size();
}

double SpineSegmentation::calculateVariance(const std::vector<double>& values, double mean) const {
	double variance = 0.0;
	for (double value : values) {
		variance += std::pow(value - mean, 2);
	}
	return variance / values.size();
}

double SpineSegmentation::calculateCorrelation(const std::vector<double>& values, int lag) const {
	std::vector<double> original(values.begin(), values.end() - lag);
	std::vector<double> lagged(values.begin() + lag, values.end());

	double meanOriginal = calculateMean(original);
	double meanLagged = calculateMean(lagged);

	double covariance = calculateCovariance(original, lagged, meanOriginal, meanLagged);
	double varianceOriginal = calculateVariance(original, meanOriginal);
	double varianceLagged = calculateVariance(lagged, meanLagged);

	return covariance / std::sqrt(varianceOriginal * varianceLagged);
}

int SpineSegmentation::getLowerValueAsPos(const std::vector<double>& values, int pBegin, int pEnd) const
{
	if (pBegin < 0 || pEnd >= values.size())
	{
		return -1;
	}

	int pos;
	for (int i = pBegin; i <= pEnd; i++)
	{
		if (i == pBegin)
		{
			pos = i;
		}
		else
		{
			if (values[i] < values[pos])
			{
				pos = i;
			}
		}
	}

	return pos;
}

int SpineSegmentation::getGreaterValueAsPos(const std::vector<double>& values, int pBegin, int pEnd) const
{
	if (pBegin < 0 || pEnd >= values.size())
	{
		return -1;
	}

	int pos;
	for (int i = pBegin; i <= pEnd; i++)
	{
		if (i == pBegin)
		{
			pos = i;
		}
		else
		{
			if (values[i] > values[pos])
			{
				pos = i;
			}
		}
	}

	return pos;
}


