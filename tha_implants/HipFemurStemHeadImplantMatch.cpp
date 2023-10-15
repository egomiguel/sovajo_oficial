#include "HipFemurStemHeadImplantMatch.hpp"
#include "ImplantsException.hpp"
#include "ImplantTools.hpp"
#include <itkVersorRigid3DTransform.h>

using namespace THA::IMPLANTS;

HipFemurStemHeadImplantMatch::HipFemurStemHeadImplantMatch()
{
    isInit = false;
}

void HipFemurStemHeadImplantMatch::init(const HipFemurStemHeadImplant& pImplantHead, const HipFemurStemImplant& pImplantStem)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_FEMUR_STEM_IMPLANT_MATCH;
    }

    this->mImplantStem = pImplantStem;
	this->mImplantHead = pImplantHead;
    getRigidTransform();
    isInit = true;
}

itk::Matrix< double, 3, 3 > HipFemurStemHeadImplantMatch::GetRotationMatrix() const
{
    itk::Matrix< double, 3, 3 > rotation;

    rotation[0][0] = rotationMatrix.at <double>(0, 0);
    rotation[0][1] = rotationMatrix.at <double>(0, 1);
    rotation[0][2] = rotationMatrix.at <double>(0, 2);

    rotation[1][0] = rotationMatrix.at <double>(1, 0);
    rotation[1][1] = rotationMatrix.at <double>(1, 1);
    rotation[1][2] = rotationMatrix.at <double>(1, 2);

    rotation[2][0] = rotationMatrix.at <double>(2, 0);
    rotation[2][1] = rotationMatrix.at <double>(2, 1);
    rotation[2][2] = rotationMatrix.at <double>(2, 2);

    return rotation;
}

itk::Vector< double, 3 > HipFemurStemHeadImplantMatch::GetTranslationMatrix() const
{
    itk::Vector< double, 3 > translation;

    translation[0] = translationMatrix.at <double>(0, 0);
    translation[1] = translationMatrix.at <double>(1, 0);
    translation[2] = translationMatrix.at <double>(2, 0);

    return translation;
}

void HipFemurStemHeadImplantMatch::getRigidTransform()
{
	Point stemNeckVector = mImplantStem.getVectorNeckToHead();

	Point stemHeadVector = mImplantHead.getVectorInfSup();

	double angle = ImplantTools::getAngleBetweenVectors(stemNeckVector, stemHeadVector);

	Point rotationAxisVector = stemHeadVector.cross(stemNeckVector);
	rotationAxisVector.normalice();

	rotationMatrix = ImplantTools::getRotateMatrix(rotationAxisVector, angle);
    translationMatrix = mImplantStem.getHeadCenter().ToMatPoint() - (rotationMatrix * mImplantHead.getInsideCenterTopPoint().ToMatPoint());
}

itk::Rigid3DTransform<>::Pointer HipFemurStemHeadImplantMatch::getStemHeadTransform() const
{
	return getITKTransform(rotationMatrix, translationMatrix);
}

itk::Rigid3DTransform<>::Pointer HipFemurStemHeadImplantMatch::getITKTransform(const cv::Mat& pRotation, const cv::Mat& pTranslation) const
{
	itk::Matrix< double, 3, 3 > rotation;
	itk::Vector< double, 3 > translate;

	rotation[0][0] = pRotation.at<double>(0, 0);
	rotation[0][1] = pRotation.at<double>(0, 1);
	rotation[0][2] = pRotation.at<double>(0, 2);

	rotation[1][0] = pRotation.at<double>(1, 0);
	rotation[1][1] = pRotation.at<double>(1, 1);
	rotation[1][2] = pRotation.at<double>(1, 2);

	rotation[2][0] = pRotation.at<double>(2, 0);
	rotation[2][1] = pRotation.at<double>(2, 1);
	rotation[2][2] = pRotation.at<double>(2, 2);

	translate[0] = pTranslation.at<double>(0, 0);
	translate[1] = pTranslation.at<double>(1, 0);
	translate[2] = pTranslation.at<double>(2, 0);

	itk::Rigid3DTransform<double>::Pointer transform = itk::VersorRigid3DTransform<double>::New();
	transform->SetMatrix(rotation);
	transform->SetOffset(translate);

	return transform;
}
