#include <opencv2/dnn.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <iostream>


std::string modelPathPB = "D:\\Algorith_Project\\DNN_MODELS\\frozen_graph_TF_example_optimize.pb";
std::string modelPathPB_TXT = "D:\\Algorith_Project\\DNN_MODELS\\frozen_graph_TF_example_optimize.pbtxt";

void segmentationDNN()
{
    std::cout << "11111111111111111111111111111111" << std::endl;
    cv::dnn::Net net;// = cv::dnn::readNet(modelPathPB, modelPathPB_TXT, "TensorFlow");

    try
    {
        net = cv::dnn::readNet(modelPathPB, modelPathPB_TXT, "TensorFlow");
    }
    catch (cv::Exception& e)
    {
        const char* err_msg = e.what();
        std::cout << "exception caught: " << err_msg << std::endl;
    }

    std::cout << "222222222222222222222222222222222222222" << std::endl;
    cv::Mat img = cv::imread("D:\\Deep_Learning\\keeshond_152.jpg", cv::IMREAD_COLOR);

    std::cout << "333333333333333333333333333333333" << std::endl;
    int down_width = 128;

    int down_height = 128;

    cv::Mat resized_down;

    cv::resize(img, resized_down, cv::Size(down_width, down_height), cv::INTER_LINEAR);
    std::cout << "4444444444444444444444444444444444444" << std::endl;
    /*net.setInput(resized_down);
    std::cout << "55555555555555555555555555555555555" << std::endl;
    cv::Mat outPut = net.forward();
    std::cout << "66666666666666666666666666666666666" << std::endl;
    cv::imshow("Segmentation", outPut);*/

}