#include "HipPelvisLinerCupImplatMatch.hpp"
#include "ImplantsException.hpp"
#include "ImplantTools.hpp"
#include <itkVersorRigid3DTransform.h>

using namespace THA::IMPLANTS;

HipPelvisLinerCupImplantMatch::HipPelvisLinerCupImplantMatch()
{
    isInit = false;
}

void HipPelvisLinerCupImplantMatch::init(const HipPelvisCupImplant& pCupImplant, const HipPelvisLinerImplant& pLinerImplant)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_PELVIS_LINER_CUT_IMPLANT_MATCH;
    }

    this->mCupImplant = pCupImplant;
    this->mLinerImplant = pLinerImplant;
	getRigidTransform();
    isInit = true;
}

Point HipPelvisLinerCupImplantMatch::getLinerCenterInCup() const
{
	return mLinerCenterInCup;
}

itk::Matrix< double, 3, 3 > HipPelvisLinerCupImplantMatch::GetRotationMatrix() const
{
	itk::Matrix< double, 3, 3 > rotation;

	rotation[0][0] = mRotationMatrix.at <double>(0, 0);
	rotation[0][1] = mRotationMatrix.at <double>(0, 1);
	rotation[0][2] = mRotationMatrix.at <double>(0, 2);

	rotation[1][0] = mRotationMatrix.at <double>(1, 0);
	rotation[1][1] = mRotationMatrix.at <double>(1, 1);
	rotation[1][2] = mRotationMatrix.at <double>(1, 2);

	rotation[2][0] = mRotationMatrix.at <double>(2, 0);
	rotation[2][1] = mRotationMatrix.at <double>(2, 1);
	rotation[2][2] = mRotationMatrix.at <double>(2, 2);

	return rotation;
}

itk::Vector< double, 3 > HipPelvisLinerCupImplantMatch::GetTranslationMatrix() const
{
	itk::Vector< double, 3 > translation;

	translation[0] = mTranslationMatrix.at <double>(0, 0);
	translation[1] = mTranslationMatrix.at <double>(1, 0);
	translation[2] = mTranslationMatrix.at <double>(2, 0);

	return translation;
}

void HipPelvisLinerCupImplantMatch::getRigidTransform()
{
	auto cupBase = mCupImplant.getBasePlane();
	cupBase.reverseByPoint(mCupImplant.getTopPoint());

	Point cupVector = cupBase.getNormalVector();

	Point linerVector = mLinerImplant.getCenterToTopVector();

	double angle = ImplantTools::getAngleBetweenVectors(cupVector, linerVector);

	Point rotationAxisVector = linerVector.cross(cupVector);
	rotationAxisVector.normalice();

	mRotationMatrix = ImplantTools::getRotateMatrix(rotationAxisVector, angle);
	mTranslationMatrix = mCupImplant.getCenterOfRotationImplant().ToMatPoint() - (mRotationMatrix * mLinerImplant.getCenterOfRotationImplant().ToMatPoint());
	
	if (mLinerImplant.getExternalRadius() > mCupImplant.getInternalRadius() && mCupImplant.getInternalRadius() > 0)
	{
		double move = mLinerImplant.getExternalRadius() - mCupImplant.getInternalRadius();
		mTranslationMatrix = mTranslationMatrix - move * cupVector.ToMatPoint();
	}

	auto newLinerCenter = mRotationMatrix * mLinerImplant.getCenterOfRotationImplant().ToMatPoint() + mTranslationMatrix;
	mLinerCenterInCup = Point(newLinerCenter);
}

itk::Rigid3DTransform<>::Pointer HipPelvisLinerCupImplantMatch::getTransform() const
{
	return getITKTransform(mRotationMatrix, mTranslationMatrix);
}

itk::Rigid3DTransform<>::Pointer HipPelvisLinerCupImplantMatch::getITKTransform(const cv::Mat& pRotation, const cv::Mat& pTranslation) const
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

