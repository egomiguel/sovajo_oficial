#include "HipFemurStemImplantMatch.hpp"
#include "ImplantsException.hpp"
#include "ImplantTools.hpp"
#include <itkVersorRigid3DTransform.h>

using namespace THA::IMPLANTS;

HipFemurStemImplantMatch::HipFemurStemImplantMatch()
{
    isInit = false;
}

void HipFemurStemImplantMatch::init(const HipPelvis& pPelvis, const HipFemurStemImplant& pImplant)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_HIP_FEMUR_STEM_IMPLANT_MATCH;
    }

    this->mPelvis = pPelvis;
    this->mImplant = pImplant;
    this->mHipCenterOfRotation = pPelvis.getHipCenterOfRotation();
    getRigidTransform();
    isInit = true;
}

itk::Matrix< double, 3, 3 > HipFemurStemImplantMatch::GetRotationMatrix() const
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

itk::Vector< double, 3 > HipFemurStemImplantMatch::GetTranslationMatrix() const
{
    itk::Vector< double, 3 > translation;

    translation[0] = translationMatrix.at <double>(0, 0);
    translation[1] = translationMatrix.at <double>(1, 0);
    translation[2] = translationMatrix.at <double>(2, 0);

    return translation;
}

void HipFemurStemImplantMatch::getRigidTransform()
{
    std::vector<cv::Point3d> implantVectors;
    std::vector<cv::Point3d> femurVectors;

    Point implantX = mImplant.getVectorLatMed();
    Point implantZ = mImplant.getVectorInfSup();
    Point implantY = implantX.cross(implantZ);
    implantY.normalice();

    implantVectors.push_back(implantX.ToCVPoint());
    implantVectors.push_back(implantZ.ToCVPoint());
    implantVectors.push_back(implantY.ToCVPoint());

	Point femurX = mPelvis.getFemurOperationSide().getVectorLatMed();  //getFemurVectorLatMed(mHipCenterOfRotation);
	Point femurZ = mPelvis.getFemurOperationSide().getCanalAxisVectorInfSup();  //getFemurVectorInfSup();
	//femurZ = femurZ + Point(0.1, 0.2, 0.3); // Just for test ******************************************************
    Point femurY = femurX.cross(femurZ);
    femurY.normalice();

	//femurZ = -femurX.cross(femurY);//***********************************************
	//femurZ.normalice();

    femurVectors.push_back(femurX.ToCVPoint());
    femurVectors.push_back(femurZ.ToCVPoint());
    femurVectors.push_back(femurY.ToCVPoint());

    cv::Mat implantMatrix = cv::Mat(implantVectors.size(), 3, CV_64F, implantVectors.data());
    cv::Mat pelvisMatrix = cv::Mat(femurVectors.size(), 3, CV_64F, femurVectors.data());

    cv::Mat inverse = (implantMatrix.t()).inv();
    rotationMatrix = (pelvisMatrix.t()) * inverse;
    translationMatrix = mHipCenterOfRotation.ToMatPoint() - (rotationMatrix * mImplant.getHeadCenter().ToMatPoint());

	//Alignment with the anatomical axis

	Point ref = mImplant.getCanalAxisRodCenter();
	cv::Mat refMat = rotationMatrix * ref.ToMatPoint() + translationMatrix;
	ref = Point(refMat);

	Line anatomical(mPelvis.getFemurOperationSide().getCanalAxisVectorInfSup(), mPelvis.getFemurOperationSide().getCanalAxisPoint());
	Point refProj = anatomical.getProjectPoint(ref);
	Point translationDiff = ref - refProj;
	translationMatrix = translationMatrix - translationDiff.ToMatPoint();
}

itk::Rigid3DTransform<>::Pointer HipFemurStemImplantMatch::getStemTransform(double pStemVersionAngleDegree) const
{
	cv::Mat neckAxisMat = rotationMatrix * mImplant.getVectorNeckToHead().ToMatPoint();
	Point neckAxis = Point(neckAxisMat);

	double angle = mPelvis.getFemurVersionRadian(neckAxis);

	double refAngle = (pStemVersionAngleDegree * PI) / 180.;

	itk::Rigid3DTransform<double>::Pointer transform = itk::VersorRigid3DTransform<double>::New();
	cv::Mat transformation;
	Point canalAxis = mPelvis.getFemurOperationSide().getCanalAxisVectorInfSup();

	if (refAngle == angle)
	{
		return getITKTransform(rotationMatrix, translationMatrix);
	}
	else if (refAngle > angle)
	{
		double temp = refAngle - angle;

		if (mPelvis.getSide() == PelvisSide::RIGHT_SIDE)
		{
			transformation = ImplantTools::getRotateMatrix(-canalAxis, temp);
		}
		else
		{
			transformation = ImplantTools::getRotateMatrix(canalAxis, temp);
		}
	}
	else
	{
		double temp = angle - refAngle;

		if (mPelvis.getSide() == PelvisSide::RIGHT_SIDE)
		{
			transformation = ImplantTools::getRotateMatrix(canalAxis, temp);
		}
		else
		{
			transformation = ImplantTools::getRotateMatrix(-canalAxis, temp);
		}
	}

	cv::Mat lastRotation = transformation * rotationMatrix;

	return getITKTransform(lastRotation, translationMatrix);
}

itk::Rigid3DTransform<>::Pointer HipFemurStemImplantMatch::getStemTransform() const
{
	return getITKTransform(rotationMatrix, translationMatrix);
}


itk::Rigid3DTransform<>::Pointer HipFemurStemImplantMatch::getITKTransform(const cv::Mat& pRotation, const cv::Mat& pTranslation) const
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
