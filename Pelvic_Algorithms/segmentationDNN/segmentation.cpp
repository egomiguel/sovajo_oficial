#include <iostream>
#include <fstream>
#include <vector>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "segmentation_dnn/segmentationDNN.hpp"

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include<opencv2/opencv.hpp>



void readImage(std::string path, RawImageType::Pointer& image)
{
    try
    {
        using ReaderType = itk::ImageFileReader< RawImageType >;
        typename ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(path);
        reader->Update();
        image = reader->GetOutput();
        std::cout << "Read image!!!" <<std::endl;
    }
    catch (const itk::ExceptionObject & ex)
    {
        std::cout << ex.GetDescription() << std::endl;
    }
}


void SaveImage(SegmentationImageType::Pointer img, std::string name, bool verbose = true)
{
    std::string fullName = name + ".nrrd";
    using WriterType = itk::ImageFileWriter<SegmentationImageType>;
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
            std::cout << "Write image!!!" << std::endl;
        }
    }
    catch (const itk::ExceptionObject & ex)
    {
        std::cout << ex.what() << std::endl;
    }
}


int main()
{
    std::string model = "D:\\Algorith_Project\\TensorFloww_CPP\\Models\\frozen_graph.pb";
    std::string imgPath = "D:\\Algorith_Project\\TensorFloww_CPP\\Images\\Raw_1.nrrd";

    std::string layerInput = "x";
    std::string layerOutput = "Identity";

    RawImageType::Pointer image;
    readImage("D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\6_Good\\Raw_New.nrrd", image);

    /*std::vector<uchar> datos = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
    cv::Mat matriz(4, 4, CV_8UC1, datos.data());
    
    std::cout << matriz << std::endl;

    cv::Mat result;
    cv::resize(matriz, result, cv::Size(8, 8), cv::INTER_LINEAR);

    std::cout << "*************************" << std::endl;
    std::cout << result << std::endl;*/

    /*
    cv::Mat convertIMG;
    img.convertTo(convertIMG, CV_32F);

    if (convertIMG.isContinuous() ==  false)
        printf("---------------------No es continuo");

    std::vector<float> vec;

    vec.assign((float*)convertIMG.data, (float*)convertIMG.data + convertIMG.total()*convertIMG.channels());


    std::cout << "img cols: " << convertIMG.cols << " img rows: " << convertIMG.rows << " img chanel: "<< convertIMG.channels() << std::endl;
    std::cout << "img lengh: " << convertIMG.cols * convertIMG.rows * convertIMG.channels() << " img pos 1000: " << convertIMG.at<float>(0, 100) << std::endl;

    std::cout << "vector lengh: " << vec.size() << " vector pos 20: " << vec[100] << std::endl;

    /*for (int i = 0; i < vec.size(); i++)
    {
        if (vec[i] != 255)
        {
            int rows = i / 900;
            int cols = i % 900;
            std::cout << "Resumen en pos: " << i << "  " << vec[i] << "  " << convertIMG.at<float>(rows, cols) << std::endl;
        }
    }*/

    /*cv::imshow("Result", img);
    cv::waitKey(0);

    cv::Mat resized_down;
    cv::resize(img, resized_down, cv::Size(128, 128), cv::INTER_LINEAR);

    cv::imshow("Result", resized_down);
    cv::waitKey(0);

    cv::Mat convertIMG;
    resized_down.convertTo(convertIMG, CV_32F);

    cv::imshow("Result", convertIMG);
    cv::waitKey(0);

    cv::Mat convertIMG2;
    convertIMG.convertTo(convertIMG2, CV_8U);

    cv::imshow("Result", convertIMG2);
    cv::waitKey(0);*/
    
    
    ModelInfo input(512, 512, 1);
    ModelInfo output(512, 512, 4);

    SegmentationDNN::getVersion();
    SegmentationDNN obj(model, layerInput, layerOutput, &input, &output);
    obj.Execute(image);

    SaveImage(obj.GetPelvis(), "Pelvis_Seg");
    SaveImage(obj.GetLegs(), "Legs_Seg");

    /*std::vector<float> aa = readFile("D:\\Algorith_Project\\TensorFloww_CPP\\Segmentation\\Build\\result_cpp.txt");
    std::vector<float> bb = readFile("D:\\Algorith_Project\\TensorFloww_CPP\\Segmentation\\Build\\result_python.txt");

    std::cout << aa.size() << "   " << bb.size() << std::endl;

    for (int i = 0; i < aa.size(); i++)
    {

        float cc = aa[i] - bb[i];

        if (abs(cc) > 0.0001)
        {
            std::cout << i << std::endl;
        }
    }

    std::cout << "Fin" << std::endl;*/

    char c;
    std::cin >> c;

	return 0;
}

