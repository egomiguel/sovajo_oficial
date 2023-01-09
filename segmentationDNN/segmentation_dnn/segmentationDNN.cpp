#include "segmentationDNN.hpp"
#include "PrivateData.hpp"
#include "tensorflow/c/c_api.h"
#include <stdio.h>
#include <functional> 
#include <algorithm>
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkVector.h"
#include "itkExtractImageFilter.h"
#include "itkImageDuplicator.h"

#include<opencv2/opencv.hpp>

#include <fstream>
#include <iostream>


SegmentationDNN::SegmentationDNN(const std::string& pModelPB, const std::string& pInputLayerName, const std::string& pOutputLayerName, ModelInfo* pInput, ModelInfo* pOutput)
{
    mData = new PrivateData();
    mInputLayerName = pInputLayerName;
    mOutputLayerName = pOutputLayerName;
    mInput = pInput;
    mOutput = pOutput;

    bool result = mData->init(pModelPB);

    if (result == true)
    {
        printf("Model Fine !!!!! \n");
    }
    else
    {
        printf("Model Wrong!!!!! \n");
    }
}

SegmentationDNN::~SegmentationDNN()
{
    mData->~PrivateData();

    mData = NULL;
    delete mData;

    mInput = NULL;
    delete mInput;

    mOutput = NULL;
    delete mOutput;
}

bool SegmentationDNN::Execute(const RawImageType::Pointer& pImg, int pPelvisLabel, int pLegsLabel)
{
    /*RawImageType::PixelType* buffer = pImg->GetBufferPointer();
    std::vector<RawPixelType> temp(buffer, buffer + 512 * 512);

    std::vector<float> vectorData;
    std::transform(temp.begin(), temp.end(), std::back_inserter(vectorData), [](RawPixelType x) { return (float)x; });

    cv::Mat A(mInput->mHeight, mInput->mWidth, CV_32FC1, vectorData.data());

    cv::Mat transpose = (A.t()) / 255.0;

    std::vector<float> tempArray;
    tempArray.assign((float*)transpose.data, (float*)transpose.data + 512 * 512);

    const std::vector<std::int64_t> input_dims = { 1, mInput->mHeight, mInput->mWidth, mInput->mChannel };
    const std::vector<std::int64_t> output_dims = { mOutput->mHeight, mOutput->mWidth, mOutput->mChannel };

    std::vector<int> prediction = GetPrediction(input_dims, output_dims, tempArray);

    std::vector<float> vectorData2;
    std::transform(prediction.begin(), prediction.end(), std::back_inserter(vectorData2), [](int x) { return (float)x / 3.0; });


    cv::Mat rawImage = cv::Mat(cv::Size(mOutput->mWidth, mOutput->mHeight), CV_32F, vectorData.data());
    cv::Mat labelImage = cv::Mat(cv::Size(mOutput->mWidth, mOutput->mHeight), CV_32F, vectorData2.data());

    cv::imshow("Raw", transpose);
    cv::imshow("labels", labelImage);
    cv::waitKey(0);*/

    auto startClock = std::chrono::system_clock::now();

    SegmentationImageType::Pointer imageIn = GetGrayScale(pImg);

    mPelvis = CloneImage(imageIn);
    mLegs = CloneImage(imageIn);

    typedef itk::ExtractImageFilter< SegmentationImageType, SegmentationImageType > FilterType;
    SegmentationImageType::RegionType inputRegion = imageIn->GetLargestPossibleRegion();
    SegmentationImageType::SizeType size = inputRegion.GetSize();

    int fullSize2D = size[0] * size[1];
    int fullSize3D = size[0] * size[1] * size[2];
    int rows = size[0];
    int cols = size[1];

    SegmentationImageType::IndexType start = inputRegion.GetIndex();
    int initZ = start[2];
    int endZ = size[2] + initZ;

    /*
    //////////////////////////////////
    int myZ = size[2];

    SegmentationImageType::PixelType* imageBuffer = imageIn->GetBufferPointer();

    const std::vector<std::int64_t> input_dims = { 1, mInput->mHeight, mInput->mWidth, mInput->mChannel };
    const std::vector<std::int64_t> output_dims = { mOutput->mHeight, mOutput->mWidth, mOutput->mChannel };

    for (int i = 0; i < fullSize3D; i += fullSize2D)
    {

        std::vector<SegmentationPixelType> temp(imageBuffer + i, imageBuffer + i + fullSize2D);
        std::vector<float> vectorData;

        std::transform(temp.begin(), temp.end(), std::back_inserter(vectorData), [](SegmentationPixelType x) { return ((float)((int)x)) / 255.0; });
        cv::Mat newImage(rows, cols, CV_32FC1, vectorData.data());
        cv::Mat downImage;

        cv::resize(newImage, downImage, cv::Size(128, 128), cv::INTER_LINEAR);

        std::vector<float> tempArray;
        tempArray.assign((float*)downImage.data, (float*)downImage.data + 128 * 128);

        std::vector<uchar> prediction = GetPrediction(input_dims, output_dims, tempArray);

        cv::Mat result(128, 128, CV_8UC1, prediction.data());

        cv::Mat upImage;
        cv::resize(result, upImage, cv::Size(rows, cols), cv::INTER_LINEAR);

        std::vector<uchar> fullImage;
        fullImage.assign(upImage.data, upImage.data + 512 * 512);

        SegmentationImageType::PixelType* pelvisBuffer = mPelvis->GetBufferPointer();
        SegmentationImageType::PixelType* legsBuffer = mLegs->GetBufferPointer();

        int zAxis = i / fullSize2D;
        for (int j = 0; j < size[1]; j++)
        {
            for (int k = 0; k < size[0]; k++)
            {
                int indexBig = zAxis * size[1] * size[0] + j * size[0] + k;
                int indexShort = j * size[0] + k;

                int var = fullImage[indexShort];

                if (var == pPelvisLabel)
                {
                    pelvisBuffer[indexBig] = 1;
                }
                else
                {
                    pelvisBuffer[indexBig] = 0;
                }

                if (var == pLegsLabel)
                {
                    legsBuffer[indexBig] = 1;
                }
                else
                {
                    legsBuffer[indexBig] = 0;
                }
            }
        }


        /*std::vector<float> prediction2;
        std::transform(prediction.begin(), prediction.end(), std::back_inserter(prediction2), [](int x) { return (float)x / 3.0; });
        cv::Mat labelImage = cv::Mat(cv::Size(mOutput->mWidth, mOutput->mHeight), CV_32F, prediction2.data());*/

       /* cv::Mat img2;
        upImage.convertTo(img2, CV_32FC1);
        img2 = img2 / 3.0;

        cv::imshow("Raw", downImage);
        cv::imshow("labels", img2);
        cv::waitKey(0);
    }





    std::cout << "Finish" << std::endl;
    ///////////////////////////////////

    return true;
    */
    size[2] = 1;

    for (int i = initZ; i < endZ; i++)
    {
        start[2] = i;
        SegmentationImageType::RegionType desiredRegion;
        desiredRegion.SetSize(size);
        desiredRegion.SetIndex(start);

        FilterType::Pointer filter = FilterType::New();
        filter->SetDirectionCollapseToSubmatrix();
        filter->SetExtractionRegion(desiredRegion);
        filter->SetInput(imageIn);
        filter->Update();

        SegmentationImageType::Pointer slice = filter->GetOutput();
        SegmentationImageType::PixelType* imageBuffer = slice->GetBufferPointer();

        std::vector<SegmentationPixelType> temp(imageBuffer, imageBuffer + fullSize2D);
        std::vector<float> vectorData;

        std::transform(temp.begin(), temp.end(), std::back_inserter(vectorData), [](SegmentationPixelType x) { return ((float)((int)x)) / 255.0; });

        ////////////////////////////////////////will be remove

        /*cv::Mat A(mInput->mHeight, mInput->mWidth, CV_32FC1, vectorData.data());
        cv::Mat transpose = A.t();

        std::vector<float> tempArray;
        tempArray.assign((float*)transpose.data, (float*)transpose.data + 512 * 512);*/

        /////////////////////////////////////////

        const std::vector<std::int64_t> input_dims = { 1, mInput->mHeight, mInput->mWidth, mInput->mChannel };
        const std::vector<std::int64_t> output_dims = { mOutput->mHeight, mOutput->mWidth, mOutput->mChannel };

        std::vector<uchar> prediction = GetPrediction(input_dims, output_dims, vectorData);

        SegmentationImageType::PixelType* pelvisBuffer = mPelvis->GetBufferPointer();
        SegmentationImageType::PixelType* legsBuffer = mLegs->GetBufferPointer();

        int zAxis = i;
        for (int j = 0; j < size[1]; j++)
        {
            for (int k = 0; k < size[0]; k++)
            {
                int indexBig = zAxis * size[1] * size[0] + j * size[0] + k;
                int indexShort = j * size[0] + k;

                int var = prediction[indexShort];

                if (var == pPelvisLabel)
                {
                    pelvisBuffer[indexBig] = 1;
                }
                else
                {
                    pelvisBuffer[indexBig] = 0;
                }

                if (var == pLegsLabel)
                {
                    legsBuffer[indexBig] = 1;
                }
                else
                {
                    legsBuffer[indexBig] = 0;
                }
            }
        }

        //cv::Mat rawImage = cv::Mat(cv::Size(mOutput->mWidth, mOutput->mHeight), CV_32F, vectorData.data());

        /////////////////////////////////////////// will be remove it
        //std::vector<float> prediction2;
        //std::transform(prediction.begin(), prediction.end(), std::back_inserter(prediction2), [](int x) { return (float)x / 3.0; });
        /////////////////////////////////////////////

        //cv::Mat labelImage = cv::Mat(cv::Size(mOutput->mWidth, mOutput->mHeight), CV_32F, prediction2.data());


        //cv::imshow("Raw", rawImage);
        //cv::imshow("labels", labelImage);
        //cv::waitKey(0);


        //return cv::Mat(cv::Size(mOutput->mWidth, mOutput->mHeight), CV_32F, finalData.data());
    }

    auto endClock = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    return true;
}

std::vector<UINT8> SegmentationDNN::GetPrediction(const std::vector<std::int64_t>& pInput_dims, const std::vector<std::int64_t>& pOutput_dims, const std::vector<float>& pData)
{
    std::vector<UINT8> resultLabels;

    TF_Output t0 = { TF_GraphOperationByName(mData->getGraph(), mInputLayerName.c_str()), 0 };
    TF_Output t1 = { TF_GraphOperationByName(mData->getGraph(), mOutputLayerName.c_str()), 0 };

    if (t0.oper == NULL) {
        printf("ERROR: Failed TF_GraphOperationByName Input\n");
    }

    if (t1.oper == NULL) {
        printf("ERROR: Failed TF_GraphOperationByName Output\n");
    }

    const std::vector<TF_Output> input_ops = { t0 };
    const std::vector<TF_Tensor*> input_tensors = { mData->CreateTensor(TF_FLOAT, pInput_dims, pData.data(), sizeof(float) * pData.size()) };

    const std::vector<TF_Output> out_ops = { t1 };
    std::vector<TF_Tensor*> output_tensors = { nullptr };

    auto code = mData->Execute(input_ops, input_tensors, out_ops, output_tensors);

    int size = pOutput_dims[0] * pOutput_dims[1];
    int channel = pOutput_dims[2];

    if (code == TF_OK)
    {
        printf("Execution Fine!!!! \n");

        auto result = mData->GetTensorData<float>(output_tensors[0]);

        float * rawData = result.data();

        for (int i = 0; i < size; i++)
        {
            int pos = i * channel;
            float maxElement = result[pos];
            int maxPos = pos;

            for (int k = 0; k < channel; k++)
            {
                if (result[pos + k] >= maxElement)
                {
                    maxElement = result[pos + k];
                    maxPos = k;
                }
            }

            resultLabels.push_back(maxPos);
        }

        printf("Return labels!!!! \n");

        return resultLabels;
    }
    else
    {
        printf("Execution Wrong!!!!! \n");
        return resultLabels;
    }
}

SegmentationImageType::Pointer SegmentationDNN::GetGrayScale(RawImageType::Pointer pImage)
{
    using RescaleType = itk::RescaleIntensityImageFilter<RawImageType, RawImageType>;
    RescaleType::Pointer rescale = RescaleType::New();
    rescale->SetInput(pImage);
    rescale->SetOutputMinimum(0);
    rescale->SetOutputMaximum(itk::NumericTraits<uint8_t>::max());

    using FilterType = itk::CastImageFilter<RawImageType, SegmentationImageType>;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(rescale->GetOutput());
    filter->Update();
    return filter->GetOutput();
}

SegmentationImageType::Pointer SegmentationDNN::CloneImage(const SegmentationImageType::Pointer input)
{
    using DuplicatorType = itk::ImageDuplicator<SegmentationImageType>;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(input);
    duplicator->Update();
    return duplicator->GetOutput();
}

SegmentationImageType::Pointer SegmentationDNN::GetPelvis() const
{
    return mPelvis;
}

SegmentationImageType::Pointer SegmentationDNN::GetLegs() const
{
    return mLegs;
}

void SegmentationDNN::getVersion()
{
    printf("Hello from TensorFlow C library version %s\n", TF_Version());
}