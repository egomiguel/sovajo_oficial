#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkGDCMImageIO.h"
#include "itkDirectory.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include <itkRigid3DTransform.h>
#include <itkVersorRigid3DTransform.h>
#include <itkMatrix.h>
#include <opencv2/calib3d/calib3d.hpp>
//#include <pcl/point_cloud.h>
//#include <pcl/point_types.h>
#include <random>
#include "types_test.hpp"
#include <fstream>
#include <set>

#include <vtkImageActor.h>
#include <vtkImageData.h>
#include <vtkImageMapper3D.h>
#include <vtkImageStencil.h>
#include <vtkImageStencilData.h>
#include <vtkImageToImageStencil.h>
#include <vtkInteractorStyleImage.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <itkBinaryFillholeImageFilter.h>
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkVector.h"
#include "itkExtractImageFilter.h"

namespace Test {

    template<typename ImageType>
    void Get2DSlice(const typename ImageType::Pointer& imageIn, typename ImageType::Pointer& imageOut, const int sliceNumber)
    {
        typedef ImageType InputImageType;
        typedef ImageType SliceImageType;
        typedef itk::ExtractImageFilter< InputImageType, SliceImageType > FilterType;
        FilterType::Pointer filter = FilterType::New();
        filter->SetDirectionCollapseToSubmatrix();
        InputImageType::RegionType inputRegion = imageIn->GetLargestPossibleRegion();
        InputImageType::SizeType size = inputRegion.GetSize();
        size[2] = 1;
        InputImageType::IndexType start = inputRegion.GetIndex();
        start[2] = sliceNumber;

        InputImageType::RegionType desiredRegion;
        desiredRegion.SetSize(size);
        desiredRegion.SetIndex(start);
        filter->SetExtractionRegion(desiredRegion);
        filter->SetInput(imageIn);
        try
        {
            filter->Update();
        }
        catch (const itk::ExceptionObject & ex)
        {
            std::cout << ex.GetDescription() << std::endl;
        }
        imageOut = filter->GetOutput();
    }

    template<typename ImageType>
    void CutVolume(const typename ImageType::Pointer& imageIn, typename ImageType::Pointer& imageOut, const int sliceNumber)
    {
        ImageType::RegionType inputRegion = imageIn->GetLargestPossibleRegion();
        ImageType::SizeType dims = inputRegion.GetSize();
        ImageType::PixelType* imageBuffer = imageIn->GetBufferPointer();

        for (int k = sliceNumber; k < dims[2]; k++)
        {
            for (int j = 0; j < dims[1]; j++)
            {
                for (int i = 0; i < dims[0]; i++)
                {
                    int index = k * dims[1] * dims[0] + j * dims[0] + i;
                    //imageBuffer[index] = value;
                }
            }
        }
    }

    template<typename ImageType>
    void ChangeSlice(typename ImageType::Pointer& pImage, int stillSlice, int value)
    {
        ImageType::RegionType inputRegion = pImage->GetLargestPossibleRegion();
        ImageType::SizeType dims = inputRegion.GetSize();
        ImageType::PixelType* imageBuffer = pImage->GetBufferPointer();

        for (int k = 0; k < stillSlice; k++)
        {
            for (int j = 0; j < dims[1]; j++)
            {
                for (int i = 0; i < dims[0]; i++)
                {
                    int index = k * dims[1] * dims[0] + j * dims[0] + i;
                    imageBuffer[index] = value;
                }
            }
        }
    }

    template<typename ImageType>
    void ChangePixel(typename ImageType::Pointer& pImage, int oldValue, int newValue)
    {
        ImageType::RegionType inputRegion = pImage->GetLargestPossibleRegion();
        ImageType::SizeType dims = inputRegion.GetSize();
        ImageType::PixelType* imageBuffer = pImage->GetBufferPointer();

        for (int k = 0; k < dims[2]; k++)
        {
            for (int j = 0; j < dims[1]; j++)
            {
                for (int i = 0; i < dims[0]; i++)
                {
                    int index = k * dims[1] * dims[0] + j * dims[0] + i;
                    int var = (int)imageBuffer[index];
                    if (var == oldValue)
                    {
                        imageBuffer[index] = newValue;
                    }
                }
            }
        }
    }

    void Generate_Full_Seg(const itk::Image<uint8_t, 3>::Pointer pSeg, itk::Image<uint8_t, 3>::Pointer& pMaso, int initSlice = 0)
    {
        constexpr unsigned int Dimension = 3;
        using ImageType = itk::Image<uint8_t, Dimension>;

        ImageType::RegionType inputRegion = pSeg->GetLargestPossibleRegion();
        ImageType::SizeType dims = inputRegion.GetSize();
        ImageType::PixelType* imageBufferMaso = pMaso->GetBufferPointer();
        ImageType::PixelType* imageBufferSeg = pSeg->GetBufferPointer();

        for (int k = initSlice; k < dims[2]; k++)
        {
            for (int j = 0; j < dims[1]; j++)
            {
                for (int i = 0; i < dims[0]; i++)
                {
                    int index = k * dims[1] * dims[0] + j * dims[0] + i;
                    int value = (int)imageBufferSeg[index];
                    if (value > 0)
                    {
                        value = value + 1;
                        if (value > 3)
                        {
                            value = 3;
                        }

                        imageBufferMaso[index] = value;
                    }
                }
            }
        }
    }

    itk::Image<uint8_t, 3>::Pointer DICOM_TO_GrayScale(itk::Image<int16_t, 3>::Pointer pImage)
    {
        constexpr unsigned int Dimension = 3;

        using InputImageType = itk::Image<int16_t, Dimension>;
        using OutputImageType = itk::Image<uint8_t, Dimension>;

        using RescaleType = itk::RescaleIntensityImageFilter<InputImageType, InputImageType>;
        RescaleType::Pointer rescale = RescaleType::New();
        rescale->SetInput(pImage);
        rescale->SetOutputMinimum(0);
        rescale->SetOutputMaximum(itk::NumericTraits<uint8_t>::max());

        using FilterType = itk::CastImageFilter<InputImageType, OutputImageType>;
        FilterType::Pointer filter = FilterType::New();
        filter->SetInput(rescale->GetOutput());
        filter->Update();
        return filter->GetOutput();
    }

    std::vector<itk::Image<uint8_t, 3>::Pointer> DICOM_TO_GrayScale_Slices_CV(itk::Image<int16_t, 3>::Pointer pImage, bool isMask, int initSlice = 0, int totalLabel = 3)
    {
        constexpr unsigned int Dimension = 3;

        using InputImageType = itk::Image<int16_t, Dimension>;
        using OutputImageType = itk::Image<uint8_t, Dimension>;
        OutputImageType::Pointer result;

        InputImageType::RegionType inputRegion = pImage->GetLargestPossibleRegion();
        InputImageType::SizeType dims = inputRegion.GetSize();

        std::set<int> values;

        if (isMask == false)
        {
            using RescaleType = itk::RescaleIntensityImageFilter<InputImageType, InputImageType>;
            RescaleType::Pointer rescale = RescaleType::New();
            rescale->SetInput(pImage);
            rescale->SetOutputMinimum(0);
            rescale->SetOutputMaximum(itk::NumericTraits<uint8_t>::max());

            using FilterType = itk::CastImageFilter<InputImageType, OutputImageType>;
            FilterType::Pointer filter = FilterType::New();
            filter->SetInput(rescale->GetOutput());
            filter->Update();
            result = filter->GetOutput();
        }
        else
        {
            using FilterType = itk::CastImageFilter<InputImageType, OutputImageType>;
            FilterType::Pointer filter = FilterType::New();
            filter->SetInput(pImage);
            filter->Update();
            result = filter->GetOutput();

            OutputImageType::PixelType* imageBuffer = result->GetBufferPointer();

            for (int k = 0; k < dims[2]; k++)
            {
                for (int j = 0; j < dims[1]; j++)
                {
                    for (int i = 0; i < dims[0]; i++)
                    {
                        int index = k * dims[1] * dims[0] + j * dims[0] + i;
                        int var = (int)imageBuffer[index];

                        if (var > totalLabel - 1)
                        {
                            var = totalLabel - 1;
                        }
                        imageBuffer[index] = var;
                        values.insert((int)imageBuffer[index]);
                    }
                }
            }

            /*
            OutputImageType::Pointer grayImage = filter->GetOutput();

            result = OutputImageType::New();

            OutputImageType::IndexType start;

            start[0] = 0;  // first index on X
            start[1] = 0;  // first index on Y
            start[2] = 0;  // first index on Z

            OutputImageType::SizeType  size;

            size[0] = dims[0];  // size along X
            size[1] = dims[1];  // size along Y
            size[2] = dims[2];  // size along Z

            OutputImageType::RegionType region;

            region.SetSize(size);
            region.SetIndex(start);

            result->SetRegions(region);
            result->Allocate();
            result->FillBuffer(0);

            OutputImageType::PixelType* grayBuffer = grayImage->GetBufferPointer();

            OutputImageType::PixelType* imageBuffer = result->GetBufferPointer();

            if (dims[0] <= 512 && dims[1] <= 512)
            {
                for (int k = 0; k < dims[2]; k++)
                {
                    for (int j = 0; j < dims[1]; j++)
                    {
                        for (int i = 0; i < dims[0]; i++)
                        {
                            int index = k * dims[1] * dims[0] + j * dims[0] + i;
                            int var = (int)grayBuffer[index];

                            if (var > 2)
                            {
                                var = 2;
                            }
                            imageBuffer[index] = var;
                            values.insert((int)imageBuffer[index]);
                        }
                    }
                }
            }
            */
        }

        int z = dims[2];
        int size = dims[0] * dims[1];
        std::vector< itk::Image<uint8_t, 3>::Pointer > allSlices;
        if (initSlice < 0)
        {
            initSlice = 0;
        }

        for (int i = initSlice; i < z; i++)
        {
            OutputImageType::Pointer imageOut;
            Get2DSlice<OutputImageType>(result, imageOut, i);
            allSlices.push_back(imageOut);

            /*OutputImageType::PixelType* buffer = imageOut->GetBufferPointer();

            cv::Mat m(dims[0], dims[1], CV_8UC1);
            memcpy(m.data, buffer, size * sizeof(uchar));
            allSlices.push_back(m);*/
        }

        for (std::set<int>::iterator it1 = values.begin(); it1 != values.end(); ++it1)
        {
            std::cout << *it1 << std::endl;
        }

        return allSlices;
    }

    void Generate_DeepLearning_csv(itk::Image<uint8_t, 3>::Pointer pImage, std::string path)
    {
        using ImageType = itk::Image<uint8_t, 3>;
        ImageType::RegionType inputRegion = pImage->GetLargestPossibleRegion();
        ImageType::SizeType dims = inputRegion.GetSize();
        ImageType::PixelType* buffer = pImage->GetBufferPointer();

        ofstream myfile(path);

        int total = dims[0] * dims[1] * dims[2];

        for (int k = 0; k < dims[2]; k++)
        {
            for (int j = 0; j < dims[1]; j++)
            {
                for (int i = 0; i < dims[0]; i++)
                {
                    int index = k * dims[1] * dims[0] + j * dims[0] + i;
                    int var = (int)buffer[index];
                    myfile << var;
                    if (k < total - 1)
                    {
                        myfile << ",";
                    }
                }
            }
        }

        myfile.close();
    }

    itk::Image<uint8_t, 3>::Pointer Generate_DeepLearning_512x1024_Mask(itk::Image<int16_t, 3>::Pointer pImage)
    {
        std::set<int> values;
        int topZ = 1024;

        using InputImageType = itk::Image<int16_t, 3>;
        using OutputImageType = itk::Image<uint8_t, 3>;

        using FilterType = itk::CastImageFilter<InputImageType, OutputImageType>;
        FilterType::Pointer filter = FilterType::New();
        filter->SetInput(pImage);
        filter->Update();
        itk::Image<uint8_t, 3>::Pointer grayImage = filter->GetOutput();

        using ImageType = itk::Image<uint8_t, 3>;

        ImageType::Pointer image = ImageType::New();

        ImageType::IndexType start;

        start[0] = 0;  // first index on X
        start[1] = 0;  // first index on Y
        start[2] = 0;  // first index on Z

        ImageType::SizeType  size;

        size[0] = 512;  // size along X
        size[1] = 512;  // size along Y
        size[2] = topZ;  // size along Z

        ImageType::RegionType region;

        region.SetSize(size);
        region.SetIndex(start);

        image->SetRegions(region);
        image->Allocate();
        image->FillBuffer(0);

        ImageType::RegionType inputRegion = grayImage->GetLargestPossibleRegion();
        ImageType::SizeType dims = inputRegion.GetSize();
        ImageType::PixelType* grayBuffer = grayImage->GetBufferPointer();

        ImageType::PixelType* imageBuffer = image->GetBufferPointer();

        if (dims[0] <= 512 && dims[1] <= 512 && dims[2] <= topZ)
        {
            for (int k = 0; k < dims[2]; k++)
            {
                for (int j = 0; j < dims[1]; j++)
                {
                    for (int i = 0; i < dims[0]; i++)
                    {
                        int index = k * dims[1] * dims[0] + j * dims[0] + i;
                        int var = (int)grayBuffer[index];
                        values.insert(var);
                        if (var > 2)
                        {
                            var = 2;
                        }
                        imageBuffer[index] = var;
                    }
                }
            }
        }
        else if (dims[2] > topZ)
        {
            int myInit = dims[2] - topZ;

            for (int k = myInit; k < dims[2]; k++)
            {
                for (int j = 0; j < dims[1]; j++)
                {
                    for (int i = 0; i < dims[0]; i++)
                    {
                        int indexGray = k * dims[1] * dims[0] + j * dims[0] + i;
                        int index = (k - myInit) * dims[1] * dims[0] + j * dims[0] + i;

                        int var = (int)grayBuffer[indexGray];
                        values.insert(var);
                        if (var > 2)
                        {
                            var = 2;
                        }
                        imageBuffer[index] = var;
                    }
                }
            }
        }

        for (std::set<int>::iterator it1 = values.begin(); it1 != values.end(); ++it1)
        {
            std::cout << *it1 << std::endl;
        }

        return image;
    }

    itk::Image<uint8_t, 3>::Pointer DICOM_TO_GrayScale_DeepLearning_512x1024(itk::Image<int16_t, 3>::Pointer pImage)
    {
        int topZ = 1024;

        itk::Image<uint8_t, 3>::Pointer grayImage = DICOM_TO_GrayScale(pImage);

        using ImageType = itk::Image<uint8_t, 3>;

        ImageType::Pointer image = ImageType::New();

        ImageType::IndexType start;

        start[0] = 0;  // first index on X
        start[1] = 0;  // first index on Y
        start[2] = 0;  // first index on Z

        ImageType::SizeType  size;

        size[0] = 512;  // size along X
        size[1] = 512;  // size along Y
        size[2] = topZ;  // size along Z

        ImageType::RegionType region;

        region.SetSize(size);
        region.SetIndex(start);

        image->SetRegions(region);
        image->Allocate();
        image->FillBuffer(0);

        ImageType::RegionType inputRegion = grayImage->GetLargestPossibleRegion();
        ImageType::SizeType dims = inputRegion.GetSize();
        ImageType::PixelType* grayBuffer = grayImage->GetBufferPointer();

        ImageType::PixelType* imageBuffer = image->GetBufferPointer();

        if (dims[0] <= 512 && dims[1] <= 512 && dims[2] <= topZ)
        {
            for (int k = 0; k < dims[2]; k++)
            {
                for (int j = 0; j < dims[1]; j++)
                {
                    for (int i = 0; i < dims[0]; i++)
                    {
                        int index = k * dims[1] * dims[0] + j * dims[0] + i;
                        imageBuffer[index] = grayBuffer[index];
                    }
                }
            }
        }
        else if (dims[2] > topZ)
        {
            int myInit = dims[2] - topZ;

            for (int k = myInit; k < dims[2]; k++)
            {
                for (int j = 0; j < dims[1]; j++)
                {
                    for (int i = 0; i < dims[0]; i++)
                    {
                        int indexGray = k * dims[1] * dims[0] + j * dims[0] + i;
                        int index = (k - myInit) * dims[1] * dims[0] + j * dims[0] + i;

                        imageBuffer[index] = grayBuffer[indexGray];
                    }
                }
            }
        }

        return image;
    }

    int ErodeImage(std::string path)
    {
        vtkNew<vtkImageToImageStencil> imageToImageStencil;
        return EXIT_SUCCESS;
    }

    void myPrint(std::string str)
    {
        std::cout << str << std::endl;
    }

    template<typename ImageType>
    bool UniformSpacing(const std::vector<std::string>& pFiles)
    {
        auto it1 = pFiles.begin();
        auto it2 = pFiles.end();
        float previous, current, distPrev = -1, distCurrent = -1, thickness = 0;
        bool sameThickness = true;
        for (; it1 != it2; ++it1)
        {
            using ReaderImage = itk::ImageFileReader<ImageType>;
            typename ReaderImage::Pointer reader = ReaderImage::New();
            reader->SetFileName(*it1);
            reader->Update();
            ImageType::Pointer image = reader->GetOutput();
            ImageType::PointType index = image->GetOrigin();

            if (it1 == pFiles.begin())
            {
                previous = index[2];
                current = index[2];
            }
            else
            {
                current = index[2];
                distCurrent = abs(current - previous);
                previous = current;

                if (distPrev < 0)
                {
                    distPrev = distCurrent;
                }

                thickness = abs(distCurrent - distPrev);
                distPrev = distCurrent;

                if (thickness > 0.05)
                {
                    sameThickness = false;
                    break;
                }
            }
        }

        return sameThickness;
    }

    template<typename ImageType>
    typename ImageType::Pointer ReadSeries(const std::string& dirName, const std::string& seriesName = "")
    {
        using NamesGeneratorType = itk::GDCMSeriesFileNames;
        NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

        nameGenerator->SetUseSeriesDetails(true);
        nameGenerator->AddSeriesRestriction("0008|0021");
        nameGenerator->SetGlobalWarningDisplay(false);
        nameGenerator->SetDirectory(dirName);

        typename ImageType::Pointer res = nullptr;

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

                bool sameThickness = UniformSpacing<ImageType>(fileNames);

                if (sameThickness == false)
                {
                    printf("Warning: Series have not uniform spacing. Possible data loss.");
                }

                using ReaderType = itk::ImageSeriesReader<ImageType>;
                typename ReaderType::Pointer reader = ReaderType::New();
                using ImageIOType = itk::GDCMImageIO;
                ImageIOType::Pointer dicomIO = ImageIOType::New();
                reader->SetImageIO(dicomIO);
                reader->SetFileNames(fileNames);
                // reader->ForceOrthogonalDirectionOff(); // properly read CTs with gantry tilt

                reader->Update();
                // Read only first series
                //typename ImageType::Pointer res = reader->GetOutput();
                res = reader->GetOutput();

                std::cout << "Load image " << std::endl;
            }
        }
        catch (itk::ExceptionObject& ex) {
            std::cout << "Can not Load image " << std::endl;
            std::cout << ex << std::endl;
            return res;

        }

        return res;
    }

    template<typename ImageType>
    void SaveImage(typename ImageType::Pointer img, std::string name, bool verbose = true)
    {
        std::string fullName = name + ".nrrd";
        using WriterType = itk::ImageFileWriter<ImageType>;
        typename WriterType::Pointer writer = WriterType::New();

        writer->SetFileName(fullName);
        writer->SetInput(img);

        if (verbose == true)
        {
            myPrint("Writing the image");
        }

        try
        {
            writer->Update();
            if (verbose == true)
            {
                myPrint("Write image!!!");
            }
        }
        catch (const itk::ExceptionObject & ex)
        {
            std::cout << ex.what() << std::endl;
        }
    }


    template<typename ImageType>
    void readImage(std::string path, typename ImageType::Pointer& image)
    {
        try
        {
            using ReaderType = itk::ImageFileReader< ImageType >;
            typename ReaderType::Pointer reader = ReaderType::New();
            reader->SetFileName(path);
            reader->Update();
            image = reader->GetOutput();
            myPrint("Read image!!!");
        }
        catch (const itk::ExceptionObject & ex)
        {
            std::cout << ex.GetDescription() << std::endl;
        }
    }

    template<typename ImageType>
    void getImageSize(typename ImageType::Pointer& image)
    {
        ImageType::SizeType  size;
        size = image->GetBufferedRegion().GetSize();
        std::cout << "Size: " << size << std::endl;
    }

    template<typename ImageType>
    typename ImageType::Pointer FillHole(typename ImageType::Pointer img)
    {
        using FillholeFilterType = itk::BinaryFillholeImageFilter<ImageType>;
        FillholeFilterType::Pointer fillHoleFilter = FillholeFilterType::New();
        fillHoleFilter->SetInput(img);
        fillHoleFilter->SetForegroundValue(1);

        try
        {
            fillHoleFilter->Update();
            myPrint("Fill!!!");
        }
        catch (const itk::ExceptionObject & ex)
        {
            std::cout << ex.what() << std::endl;
        }
        return fillHoleFilter->GetOutput();
    }

    std::vector<PointTypeITK> readITKPointFile(std::string point_path)
    {
        std::vector<PointTypeITK> points;
        std::ifstream infile(point_path);
        double a, b, c;

        while (infile >> a >> b >> c)
        {
            PointTypeITK itkPoint;
            itkPoint[0] = a;
            itkPoint[1] = b;
            itkPoint[2] = c;
            points.push_back(itkPoint);
        }
        infile.close();
        std::cout << "Points ITK Size: " << points.size() << std::endl;
        return points;
    }

    std::vector<PointTypeITK> readITKPointFile(std::vector<std::string> point_path)
    {
        std::vector<PointTypeITK> points;
        for (int i = 0; i < point_path.size(); i++)
        {
            std::ifstream infile(point_path[i]);
            double a, b, c;

            while (infile >> a >> b >> c)
            {
                PointTypeITK itkPoint;
                itkPoint[0] = a;
                itkPoint[1] = b;
                itkPoint[2] = c;
                points.push_back(itkPoint);
            }
            infile.close();
        }
        std::cout << "Points ITK Size: " << points.size() << std::endl;
        return points;
    }

    cv::Mat Rigid3DTransformToCV(itk::Rigid3DTransform<>::Pointer transform)
    {
        itk::Matrix< double, 3, 3 > rotation = transform->GetMatrix();
        itk::Vector< double, 3 > translate = transform->GetOffset();

        double * matrix = new double[16];
        matrix[0] = rotation[0][0];
        matrix[1] = rotation[0][1];
        matrix[2] = rotation[0][2];
        matrix[3] = translate[0];

        matrix[4] = rotation[1][0];
        matrix[5] = rotation[1][1];
        matrix[6] = rotation[1][2];
        matrix[7] = translate[1];

        matrix[8] = rotation[2][0];
        matrix[9] = rotation[2][1];
        matrix[10] = rotation[2][2];
        matrix[11] = translate[2];

        matrix[12] = 0;
        matrix[13] = 0;
        matrix[14] = 0;
        matrix[15] = 1;

        cv::Mat matrixCV = cv::Mat(4, 4, CV_64FC1, matrix);

        return matrixCV;
    }

    cv::Mat Rigid3DTransformToCVRotation(itk::Rigid3DTransform<>::Pointer transform)
    {
        itk::Matrix< double, 3, 3 > rotation = transform->GetMatrix();

        double * matrix = new double[9];
        matrix[0] = rotation[0][0];
        matrix[1] = rotation[0][1];
        matrix[2] = rotation[0][2];

        matrix[3] = rotation[1][0];
        matrix[4] = rotation[1][1];
        matrix[5] = rotation[1][2];

        matrix[6] = rotation[2][0];
        matrix[7] = rotation[2][1];
        matrix[8] = rotation[2][2];

        cv::Mat matrixCV = cv::Mat(3, 3, CV_64FC1, matrix);

        return matrixCV;
    }

    cv::Mat Rigid3DTransformToCVTranslation(itk::Rigid3DTransform<>::Pointer transform)
    {
        itk::Vector< double, 3 > translate = transform->GetOffset();

        double * matrix = new double[3];

        matrix[0] = translate[0];
        matrix[1] = translate[1];
        matrix[2] = translate[2];

        cv::Mat matrixCV = cv::Mat(3, 1, CV_64FC1, matrix);

        return matrixCV;
    }

    /*
    pcl::PointCloud<pcl::PointXYZ>::Ptr itkVectorToPCL(std::vector<PointTypeITK> pointsITK)
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr Points(new pcl::PointCloud<pcl::PointXYZ>);
        for (int i = 0; i < pointsITK.size(); i++)
        {
            pcl::PointXYZ temp(pointsITK[i][0], pointsITK[i][1], pointsITK[i][2]);
            Points->points.push_back(temp);
        }

        return Points;
    }
    */

    itk::Rigid3DTransform<double>::Pointer generateRigid3DTransform()
    {
        itk::Matrix< double, 3, 3 > rotation;
        itk::Vector< double, 3 > translate;
        rotation[0][0] = 1;
        rotation[0][1] = 0;
        rotation[0][2] = 0;

        rotation[1][0] = 0;
        rotation[1][1] = 1;
        rotation[1][2] = 0;

        rotation[2][0] = 0;
        rotation[2][1] = 0;
        rotation[2][2] = 1;

        translate[0] = 0;
        translate[1] = 0;
        translate[2] = 0;

        itk::Rigid3DTransform<double>::Pointer transform = itk::VersorRigid3DTransform<double>::New();
        transform->SetMatrix(rotation);
        transform->SetOffset(translate);
        return transform;
    }

    class LitePlane
    {
    public:
        cv::Point3d normalVector;
        double bias;
        cv::Point3d mPoint;
        std::string name = "default.txt";
        LitePlane(const LitePlane& X)
        {
            normalVector = X.normalVector;
            mPoint = X.mPoint;
            name = X.name;
            bias = X.bias;
        }
        LitePlane(cv::Point3d nVector, double pBias)
        {
            normalVector = nVector;
            bias = pBias;

            double squareNorm = sqrt(normalVector.dot(normalVector));
            normalVector = (normalVector / squareNorm);
            bias = (bias / squareNorm);
            double z = -bias / nVector.z;
            mPoint = cv::Point3d(0.0, 0.0, z);
        }
        void setName(std::string pName)
        {
            name = pName + ".txt";
        }

        LitePlane(cv::Point3d nVector, cv::Point3d mPoint)
        {
            normalVector = nVector;
            bias = (-1.0)*normalVector.dot(mPoint);

            double squareNorm = sqrt(normalVector.dot(normalVector));
            normalVector = (normalVector / squareNorm);
            bias = (bias / squareNorm);
            this->mPoint = mPoint;
        }

        bool isPointNear(cv::Point3d& pPoint, double distance = 0.1)
        {
            double evaluation = (normalVector.dot(pPoint) + bias) / sqrt(normalVector.dot(normalVector));
            if (fabs(evaluation) <= distance)
            {
                return true;
            }
            return false;
        }

        double eval(cv::Point3d a)
        {
            return normalVector.dot(a) + bias;
        }

        static cv::Point3d pointOnLineProj(cv::Point3d line1, cv::Point3d line2, cv::Point3d proj)
        {
            cv::Point3d directVector = line1 - line2;
            cv::Point3d diff = proj - line1;
            cv::Point3d projectOnDirector = ((diff.dot(directVector)) / (directVector.dot(directVector))) * directVector;
            cv::Point3d projection = line1 + projectOnDirector;
            return projection;
        }

        cv::Point3d getProjectionPoint(cv::Point3d pPoint)
        {
            cv::Point3d diff = pPoint - mPoint;
            cv::Point3d projectOnNormal = ((diff.dot(normalVector)) / (normalVector.dot(normalVector))) * normalVector;
            cv::Point3d projection = pPoint - projectOnNormal;
            return projection;
        }

        static double getAngleBetweenVectors(cv::Point3d a, cv::Point3d b)
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

        static double getDistanceBetweenPoints(cv::Point3d& a, cv::Point3d& b, bool square = false)
        {
            cv::Point3d diff = a - b;
            if (square == false)
            {
                return sqrt(diff.dot(diff));
            }
            else
            {
                return diff.dot(diff);
            }
        }

        static cv::Mat rotationMatrix(cv::Point3d axis, double angle)
        {
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

            cv::Mat I = cv::Mat::eye(3, 3, CV_64F);
            cv::Mat result = I + sin(angle)*rotationMatrix + (1.0 - cos(angle))*(rotationMatrix * rotationMatrix);
            return result;
        }

        static cv::Mat pointToMat(cv::Point3d pPoint)
        {
            cv::Mat normalPlaneMat(3, 1, CV_64F);
            normalPlaneMat.at <double>(0, 0) = pPoint.x;
            normalPlaneMat.at <double>(1, 0) = pPoint.y;
            normalPlaneMat.at <double>(2, 0) = pPoint.z;
            return normalPlaneMat;
        }

        static cv::Mat getRotateMatrix(cv::Point3d normalPlane)
        {
            cv::Mat normalPlaneMat(3, 1, CV_64F);
            normalPlaneMat.at <double>(0, 0) = normalPlane.x;
            normalPlaneMat.at <double>(1, 0) = normalPlane.y;
            normalPlaneMat.at <double>(2, 0) = normalPlane.z;

            cv::Point3d normalXY(0.0, 0.0, 1.0);
            cv::Point3d rotationAxis = normalXY.cross(normalPlane);
            rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
            double rotationAngle = getAngleBetweenVectors(normalXY, normalPlane);
            cv::Mat rotation_1 = rotationMatrix(rotationAxis, -rotationAngle);
            cv::Mat rotation_2 = rotationMatrix(rotationAxis, rotationAngle);

            cv::Mat rotateVector_1 = rotation_1 * normalPlaneMat;
            cv::Mat rotateVector_2 = rotation_2 * normalPlaneMat;

            cv::Mat rotate;

            double distance_1 = getDistanceBetweenPoints(cv::Point3d(rotateVector_1), normalXY);
            double distance_2 = getDistanceBetweenPoints(cv::Point3d(rotateVector_2), normalXY);

            if (distance_1 < distance_2)
            {
                rotate = rotation_1;
            }
            else
            {
                rotate = rotation_2;
            }
            return rotate;
        }

    };

    cv::Mat Rx(double angle)
    {
        cv::Mat matrix = cv::Mat::zeros(3, 3, CV_64F);

        matrix.at<double>(0, 0) = 1;
        matrix.at<double>(0, 1) = 0;
        matrix.at<double>(0, 2) = 0;

        matrix.at<double>(1, 0) = 0;
        matrix.at<double>(1, 1) = cos(angle);
        matrix.at<double>(1, 2) = -sin(angle);

        matrix.at<double>(2, 0) = 0;
        matrix.at<double>(2, 1) = sin(angle);
        matrix.at<double>(2, 2) = cos(angle);
        return matrix;
    }

    cv::Mat Ry(double angle)
    {
        cv::Mat matrix = cv::Mat::zeros(3, 3, CV_64F);

        matrix.at<double>(0, 0) = cos(angle);
        matrix.at<double>(0, 1) = 0;
        matrix.at<double>(0, 2) = sin(angle);

        matrix.at<double>(1, 0) = 0;
        matrix.at<double>(1, 1) = 1;
        matrix.at<double>(1, 2) = 0;

        matrix.at<double>(2, 0) = -sin(angle);
        matrix.at<double>(2, 1) = 0;
        matrix.at<double>(2, 2) = cos(angle);
        return matrix;
    }

    cv::Mat Rz(double angle)
    {
        cv::Mat matrix = cv::Mat::zeros(3, 3, CV_64F);

        matrix.at<double>(0, 0) = cos(angle);
        matrix.at<double>(0, 1) = -sin(angle);
        matrix.at<double>(0, 2) = 0;

        matrix.at<double>(1, 0) = sin(angle);
        matrix.at<double>(1, 1) = cos(angle);
        matrix.at<double>(1, 2) = 0;

        matrix.at<double>(2, 0) = 0;
        matrix.at<double>(2, 1) = 0;
        matrix.at<double>(2, 2) = 1;
        return matrix;
    }

    double GenerateRandomNumber(double a, double b)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> distr(a, b);
        return distr(gen);
    }

    cv::Mat GenerateRotationMatrix()
    {
        double a = 3.14 - GenerateRandomNumber(0, 6.2);
        double b = 3.14 - GenerateRandomNumber(0, 6.2);
        double c = 3.14 - GenerateRandomNumber(0, 6.2);
        cv::Mat Rotation = Rx(a) * Ry(b) * Rz(c);
        return Rotation;
    }

    cv::Mat GenerateRotationMatrix(double a, double b, double c)
    {
        a = (a * acos(-1)) / 180.0;
        b = (b * acos(-1)) / 180.0;
        c = (c * acos(-1)) / 180.0;

        cv::Mat Rotation = Rx(a) * Ry(b) * Rz(c);
        return Rotation;
    }

    cv::Mat GenerateTranslationMatrix(double a, double b, double c)
    {
        cv::Mat Translation(3, 1, CV_64F);
        Translation.at<double>(0, 0) = a;
        Translation.at<double>(1, 0) = b;
        Translation.at<double>(2, 0) = c;
        return Translation;
    }

    cv::Mat AddNoiseTranslation(cv::Mat matrix)
    {
        double a = 0.5 - GenerateRandomNumber(0, 1.0);
        double b = 0.5 - GenerateRandomNumber(0, 1.0);
        double c = 0.5 - GenerateRandomNumber(0, 1.0);
        cv::Mat Translation(3, 1, CV_64F);
        Translation.at<double>(0, 0) = matrix.at<double>(0, 0) + a;
        Translation.at<double>(1, 0) = matrix.at<double>(1, 0) + b;
        Translation.at<double>(2, 0) = matrix.at<double>(2, 0) + c;
        return Translation;
    }

    cv::Mat GenerateTranslationMatrix()
    {
        double a = 1000.0 - GenerateRandomNumber(0, 2000.0);
        double b = 1000.0 - GenerateRandomNumber(0, 2000.0);
        double c = 1000.0 - GenerateRandomNumber(0, 2000.0);
        cv::Mat Translation(3, 1, CV_64F);
        Translation.at<double>(0, 0) = a;
        Translation.at<double>(1, 0) = b;
        Translation.at<double>(2, 0) = c;
        return Translation;
    }

    itk::Point<double, 3> MatToITK(cv::Mat PointMat)
    {
        itk::Point<double, 3> myPoint;
        myPoint[0] = PointMat.at<double>(0, 0);
        myPoint[1] = PointMat.at<double>(1, 0);
        myPoint[2] = PointMat.at<double>(2, 0);
        return myPoint;
    }


    std::vector<itk::Point<double, 3>> TransformPointVector(std::vector<cv::Point3d> points, cv::Mat rotation, cv::Mat translation, bool noise = false, double scale = 1.0)
    {
        std::vector<itk::Point<double, 3>> result;

        for (int i = 0; i < points.size(); i++)
        {
            cv::Mat pointMat(3, 1, CV_64F);
            pointMat.at<double>(0, 0) = points[i].x;
            pointMat.at<double>(1, 0) = points[i].y;
            pointMat.at<double>(2, 0) = points[i].z;

            cv::Mat temp = scale * (rotation * pointMat + translation);

            if (noise == true)
            {
                double a = 0.5 - GenerateRandomNumber(0, 1.0);
                double b = 0.5 - GenerateRandomNumber(0, 1.0);
                double c = 0.5 - GenerateRandomNumber(0, 1.0);
                temp.at<double>(0, 0) = temp.at<double>(0, 0) + a;
                temp.at<double>(1, 0) = temp.at<double>(1, 0) + b;
                temp.at<double>(2, 0) = temp.at<double>(2, 0) + c;
            }
            result.push_back(MatToITK(temp));
        }

        return result;
    }

    itk::Point<double, 3> TransformOnePoint(cv::Point3d points, cv::Mat rotation, cv::Mat translation, bool noise = false, double scale = 1.0)
    {

        cv::Mat pointMat(3, 1, CV_64F);
        pointMat.at<double>(0, 0) = points.x;
        pointMat.at<double>(1, 0) = points.y;
        pointMat.at<double>(2, 0) = points.z;

        cv::Mat temp = scale * (rotation * pointMat + translation);

        if (noise == true)
        {
            double a = 5.0 - GenerateRandomNumber(0, 10.0);
            double b = 5.0 - GenerateRandomNumber(0, 10.0);
            double c = 5.0 - GenerateRandomNumber(0, 10.0);
            temp.at<double>(0, 0) = temp.at<double>(0, 0) + a;
            temp.at<double>(1, 0) = temp.at<double>(1, 0) + b;
            temp.at<double>(2, 0) = temp.at<double>(2, 0) + c;
        }

        return MatToITK(temp);
    }

    cv::Mat ITKToMat(itk::Point<double, 3> point)
    {
        cv::Mat PointMat(3, 1, CV_64F);
        PointMat.at<double>(0, 0) = point[0];
        PointMat.at<double>(1, 0) = point[1];
        PointMat.at<double>(2, 0) = point[2];
        return PointMat;
    }

    cv::Point3d ITKToCV(itk::Point<double, 3> point)
    {
        cv::Point3d myPoint(point[0], point[1], point[2]);
        return myPoint;
    }

    itk::Point<double, 3> AddNoiseITK(itk::Point<double, 3> pointITK, double maxValue = 4.0)
    {
        if (maxValue <= 0)
        {
            return pointITK;
        }

        double a = (maxValue / 2.0) - GenerateRandomNumber(0, maxValue);
        double b = (maxValue / 2.0) - GenerateRandomNumber(0, maxValue);
        double c = (maxValue / 2.0) - GenerateRandomNumber(0, maxValue);

        itk::Point<double, 3> myPoint;
        myPoint[0] = pointITK[0] + a;
        myPoint[1] = pointITK[1] + b;
        myPoint[2] = pointITK[2] + c;
        return myPoint;
    }

    void saveFileITK(const std::vector<itk::Point<double, 3>>& pPoints, std::string name)
    {
        std::string fullName = name + ".txt";
        ofstream MyFile(fullName);

        for (int i = 0; i < pPoints.size(); i++)
        {
            MyFile << pPoints[i][0] << ", " << pPoints[i][1] << ", " << pPoints[i][2] << "\n";
        }

        MyFile.close();
        std::cout << "Write file" << std::endl;
    }

    std::pair<cv::Point3d, double> fitSphere(const std::vector<cv::Point3d>& pPoints)
    {
        std::vector<cv::Point3d>::const_iterator it1, it2;
        it1 = pPoints.begin();
        it2 = pPoints.end();
        std::vector<cv::Mat> points;
        std::vector<double> bias;

        cv::Mat A(pPoints.size(), 4, CV_64F);
        int cont = 0;

        for (; it1 != it2; ++it1)
        {
            cv::Point3d temp = (*it1);
            A.at<double>(cont, 0) = temp.x;
            A.at<double>(cont, 1) = temp.y;
            A.at<double>(cont, 2) = temp.z;
            A.at<double>(cont, 3) = 1.0;
            bias.push_back(temp.dot(temp));
            cont++;
        }

        cv::Mat B = cv::Mat(bias.size(), 1, CV_64F, bias.data());

        cv::Mat dx;

        bool result = true;

        try
        {
            result = cv::solve(A, B, dx, cv::DECOMP_SVD);
        }
        catch (const std::exception& e)
        {
            std::cout << e.what() << std::endl;
        }

        if (result == false)
        {
            return std::make_pair(cv::Point3d(0, 0, 0), -1);
        }

        cv::Point3d center = (cv::Point3d(dx.at<double>(0, 0), dx.at<double>(1, 0), dx.at<double>(2, 0))) / 2.0;

        double radius = sqrt(center.dot(center) + dx.at<double>(3, 0));

        return std::make_pair(center, radius);
    }


    /*
bool compareInterval(cv::Point2d p1, cv::Point2d p2)
{
    double angle_p1 = atan2(p1.y, p1.x) * 180.0 / acos(-1);
    double angle_p2 = atan2(p2.y, p2.x) * 180.0 / acos(-1);
    if (angle_p1 < 0)
    {
        angle_p1 = 360.0 + angle_p1;
    }
    if (angle_p2 < 0)
    {
        angle_p2 = 360.0 + angle_p2;
    }
    return (angle_p1 < angle_p2);
}

double getAngle(cv::Point2d p)
{
    double angle = atan2(p.y, p.x) * 180.0 / acos(-1);
    if (angle < 0)
    {
        angle = 360.0 + angle;
    }
    return angle;
}

void angles()
{
    double X, Y;
    double rad, angle;
    std::vector<cv::Point2d> points;
    std::vector<cv::Point3d> points3d;
    for (double i = 0.0; i < 360.0; i += 45.0)
    {
        rad = i * acos(-1) / 180.0;

        X = cos(rad);
        Y = sin(rad);
        points.push_back(cv::Point2d(X, Y));
        angle = atan2(Y, X) * 180.0 / acos(-1);
        if (angle < 0)
        {
            angle = 360.0 + angle;
        }
        std::cout << "I " << i << ": " << angle << std::endl;
    }
    std::reverse(points.begin(), points.end());
    for (int i = 0; i < points.size(); i++)
    {
        std::cout <<"Angle: " << getAngle(points[i]) << " Point: "<< points[i] << std::endl;
        points3d.push_back(cv::Point3d(points[i].x, points[i].y, 1.0));
    }
    cv::Mat rotation = cv::Mat::eye(3, 3, CV_64F);
    sort(points3d.begin(), points3d.end(), PointsComp(rotation));

    for (int i = 0; i < points3d.size(); i++)
    {
        std::cout << "Angle: " << getAngle(cv::Point2d(points3d[i].x, points3d[i].y)) << " Point: " << points3d[i] << std::endl;
    }
}
*/
/*
void savePoints(std::vector<PointTypeITK> result, std::string name)
{
    std::ofstream myfile;
    name = name + ".txt";
    myfile.open(name);

    for (int i = 0; i < result.size(); i++)
    {
        myfile << result[i] << "\n";
    }
    myfile.close();
}


pcl::PointXYZ itkToPclPoint(PointTypeITK itkPoint, pcl::PointXYZ centroid)
{
    pcl::PointXYZ P1;
    P1.x = itkPoint[0] - centroid.x;
    P1.y = itkPoint[1] - centroid.y;
    P1.z = itkPoint[2] - centroid.z;
    return P1;
}





cv::Mat PointITKToMat4(PointTypeITK P)
{
    cv::Mat matrixCV(4, 1, CV_64FC1);
    matrixCV.at<double>(0, 0) = P[0];
    matrixCV.at<double>(1, 0) = P[1];
    matrixCV.at<double>(2, 0) = P[2];
    matrixCV.at<double>(3, 0) = 1;
    return matrixCV;
}



int generateRandom(int radius)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, radius);
    return distr(gen);
}

struct PointsComp
{
    cv::Mat rotation;
    PointsComp(const cv::Mat& pRotation) : rotation(pRotation) {}

    bool operator()(cv::Point3d point1, cv::Point3d point2) const
    {
        cv::Mat pointMat1(3, 1, CV_64F);
        cv::Mat pointMat2(3, 1, CV_64F);

        pointMat1.at <double>(0, 0) = point1.x;
        pointMat1.at <double>(1, 0) = point1.y;
        pointMat1.at <double>(2, 0) = point1.z;

        pointMat2.at <double>(0, 0) = point2.x;
        pointMat2.at <double>(1, 0) = point2.y;
        pointMat2.at <double>(2, 0) = point2.z;

        cv::Mat rotateP1 = rotation * pointMat1;
        cv::Mat rotateP2 = rotation * pointMat2;

        double angle_p1 = atan2(rotateP1.at<double>(1, 0), rotateP1.at<double>(0, 0)) * 180.0 / acos(-1);
        double angle_p2 = atan2(rotateP2.at<double>(1, 0), rotateP2.at<double>(0, 0)) * 180.0 / acos(-1);
        if (angle_p1 < 0)
        {
            angle_p1 = 360.0 + angle_p1;
        }
        if (angle_p2 < 0)
        {
            angle_p2 = 360.0 + angle_p2;
        }
        return (angle_p1 < angle_p2);
    }
};
*/

/*
double executeRegistrationTest()
{
    return 0.0;

    /*
        For the registration first step process you must have 4 ct points of the femur and 4 ct points of the tibia.
        A 3d volume of the tibial knee and a 3d volume of the femural knee.

        You must create the objects FemurRegistration for femur registration process and TibiaRegistration for
        tibia registration process as it shown below.

        For the second step you must obtain 40 points of the bone to be registered.
        These points are collected by the surgeon and store in vector like this (std::vector<PointTypeITK> realPoints)

        Finally, you just need to use the method MakeRegistration(). This method receives the 40 points collected
        from the bone in std::vector<PointTypeITK>,  4 specific points of the bone (These points correspond to the 4
        points of the CT image previously used to create the object.) and the regulation parameter alpha
        that must belong to the range [1.0; 2.0], by default is 2.0.

    */
    /*
    PointTypeITK hipCenter, kneeCenter, lateralEpi, medialEpi, cortex, tibiaKnee, tubercle, lateralAnkle, medialAnkle;
    PointTypeITK hipCenter2, kneeCenter2, lateralEpi2, medialEpi2, cortex2, tibiaKnee2, tubercle2, lateralAnkle2, medialAnkle2;
    //Registration::generate3D("C:\\DICOM\\femur");
    //return;


    hipCenter[0] = 22.5;
    hipCenter[1] = -50.5;
    hipCenter[2] = 874.0;

    kneeCenter[0] = 24.0;
    kneeCenter[1] = -64.5;
    kneeCenter[2] = 769.0;

    lateralEpi[0] = -9.4;
    lateralEpi[1] = -53.6;
    lateralEpi[2] = 778.6;

    medialEpi[0] = 65.3;
    medialEpi[1] = -50.0;
    medialEpi[2] = 771.1;

    cortex[0] = 24.5;
    cortex[1] = -76.0;
    cortex[2] = 801.1;

    //////////////////////////////////
    pcl::PointXYZ centroid(26.8179, -50.95, 794.896);
    pcl::PointXYZ pCortex, pLateral, pMedial, pKnee;

    pCortex = itkToPclPoint(cortex, centroid);
    pLateral = itkToPclPoint(lateralEpi, centroid);
    pMedial = itkToPclPoint(medialEpi, centroid);
    pKnee = itkToPclPoint(kneeCenter, centroid);

    /////////////////////////////////

    std::vector<RegistrationImageType::Pointer> imgPoints;
    RegistrationImageType::Pointer knee, hip;
    readImage<RegistrationImageType>(dicom_femur, knee);
    //readImage<RegistrationImageType>(hip3d, hip);
    //imgPoints.push_back(hip);
    imgPoints.push_back(knee);

    FemurRegistration * regis = new FemurRegistration(imgPoints, hipCenter, cortex, kneeCenter, lateralEpi, medialEpi);

    pcl::PointCloud<pcl::PointXYZ>::Ptr bonePoints;
    pcl::PointCloud<pcl::PointXYZ>::Ptr redPoints(new pcl::PointCloud<pcl::PointXYZ>);
    //redPoints->points.push_back(pCortex);
    //redPoints->points.push_back(pLateral);
    //redPoints->points.push_back(pMedial);
    //redPoints->points.push_back(pKnee);

    bonePoints = regis->getPointsCT();
    //Registration::drawCloud(bonePoints, redPoints);
    //return 0.0;

    /////////////////////////////////// Generating bone points for first registration step. This code is for testing.////
    regis->setTransformMatrix();

    hipCenter2 = regis->TransformPointTest(hipCenter, true);
    cortex2 = regis->TransformPointTest(cortex, true);
    kneeCenter2 = regis->TransformPointTest(kneeCenter, true);
    lateralEpi2 = regis->TransformPointTest(lateralEpi, true);
    medialEpi2 = regis->TransformPointTest(medialEpi, true);

    //////////////////////////////////////////////////////////////////////////////////////////
    std::vector<PointTypeITK> bones_points;

    //std::cout << "Hip Center 2: " << hipCenter2 << std::endl;

    bool result = regis->MakeRegistration(bones_points, hipCenter2, cortex2, kneeCenter2, lateralEpi2, medialEpi2, 0.2);
    //bool result = regis->MakeRegistration(bones_points, tibiaKnee2, tubercle2, lateralAnkle2, medialAnkle2, 0.2);

    std::cout << "****** Result: " << result << " Error: " << regis->getError() << std::endl;
    /*
    pcl::PointCloud<pcl::PointXYZ>::Ptr tibiaPoints(new pcl::PointCloud<pcl::PointXYZ>);

    for (int i = 0; i < bones_points.size(); i++)
    {
        tibiaPoints->points.push_back(pcl::PointXYZ(bones_points[i][0], bones_points[i][1], bones_points[i][2]));
    }
    redPoints->points.push_back(pcl::PointXYZ(tibiaKnee2[0], tibiaKnee2[1], tibiaKnee2[2]));
    redPoints->points.push_back(pcl::PointXYZ(tubercle2[0], tubercle2[1], tubercle2[2]));
    //redPoints->points.push_back(pcl::PointXYZ(lateralAnkle2[0], lateralAnkle2[1], lateralAnkle2[2]));
    //redPoints->points.push_back(pcl::PointXYZ(medialAnkle2[0], medialAnkle2[1], medialAnkle2[2]));

    Registration::drawCloud(tibiaPoints, redPoints, true);

    return regis->getError();
} */


}