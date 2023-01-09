#include <string>
#include "itkImage.h"

using SegmentationPixelType = uint8_t;
using RawPixelType = int16_t;
using RawImageType = itk::Image<RawPixelType, 3>;
using SegmentationImageType = itk::Image<SegmentationPixelType, 3>;


class PrivateData;

struct ModelInfo
{
    int mWidth, mHeight, mChannel;
    ModelInfo(int pWidth, int pHeight, int pChannel) :
        mWidth(pWidth), mHeight(pHeight), mChannel(pChannel) {}
};

class SegmentationDNN
{
public:
    SegmentationDNN(const std::string& pModelPB, const std::string& pInputLayerName, const std::string& pOutputLayerName, ModelInfo* pInput, ModelInfo* pOutput);
    
    ~SegmentationDNN();

    bool Execute(const RawImageType::Pointer& pImg, int pPelvisLabel = 2, int pLegsLabel = 3);

    SegmentationImageType::Pointer GetPelvis() const;

    SegmentationImageType::Pointer GetLegs() const;

    static void getVersion();

private:
    PrivateData* mData;
    std::string mInputLayerName, mOutputLayerName;
    ModelInfo* mInput;
    ModelInfo* mOutput;
    SegmentationImageType::Pointer mPelvis, mLegs;
    SegmentationImageType::Pointer GetGrayScale(RawImageType::Pointer pImage);
    std::vector<UINT8> GetPrediction(const std::vector<std::int64_t>& pInput_dims, const std::vector<std::int64_t>& pOutput_dims, const std::vector<float>& pData);
    SegmentationImageType::Pointer CloneImage(const SegmentationImageType::Pointer input);
};
