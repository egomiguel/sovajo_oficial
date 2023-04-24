#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include "segmentation_dnn/segmentationDNN.hpp"


int main()
{
    std::string model = "D:\\Algorith_Project\\TensorFloww_CPP\\Segmentation_Body\\Source\\segmentation_dnn\\model\\pix_model.pb";

    auto startClock = std::chrono::system_clock::now();

    cv::Mat img = cv::imread("D:\\Deep_Learning\\data_human\\hombre.png", cv::IMREAD_COLOR);

    SegmentationDNN obj(model);
    obj.Execute(img);

    cv::Mat mask = obj.getBodyMask();
    cv::Mat parts = obj.getBodyParts();

    auto endClock = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    cv::Mat viewImg;
    cv::hconcat(mask, parts, viewImg);

    cv::imshow("Images", viewImg);
    cv::waitKey(0);

	return 0;
}

