#include "HipFemurImplantsMatchInfo.hpp"
#include "ImplantsException.hpp"
#include "ImplantTools.hpp"
#include <itkVersorRigid3DTransform.h>

using namespace THA::IMPLANTS;

HipFemurImplantsMatchInfo::HipFemurImplantsMatchInfo(const HipPelvis& pPelvis, const HipFemurStemImplant& pImplantStem, const itk::Rigid3DTransform<>::Pointer pImplantToBoneStemTransform)
{
	this->mPelvis = pPelvis;
	this->mImplantStem = pImplantStem;
	setStemTransform(pImplantToBoneStemTransform);
}

HipFemurImplantsMatchInfo::~HipFemurImplantsMatchInfo()
{
}

void HipFemurImplantsMatchInfo::setStemTransform(const itk::Rigid3DTransform<>::Pointer pImplantToBoneStemTransform)
{
	mRotationStem = Rigid3DTransformToCVRotation(pImplantToBoneStemTransform);
	mTranslationStem = Rigid3DTransformToCVTranslation(pImplantToBoneStemTransform);
}

itk::Rigid3DTransform<>::Pointer HipFemurImplantsMatchInfo::getITKStemTransform() const
{
	itk::Matrix< double, 3, 3 > rotation;
	itk::Vector< double, 3 > translate;

	rotation[0][0] = mRotationStem.at<double>(0, 0);
	rotation[0][1] = mRotationStem.at<double>(0, 1);
	rotation[0][2] = mRotationStem.at<double>(0, 2);

	rotation[1][0] = mRotationStem.at<double>(1, 0);
	rotation[1][1] = mRotationStem.at<double>(1, 1);
	rotation[1][2] = mRotationStem.at<double>(1, 2);

	rotation[2][0] = mRotationStem.at<double>(2, 0);
	rotation[2][1] = mRotationStem.at<double>(2, 1);
	rotation[2][2] = mRotationStem.at<double>(2, 2);

	translate[0] = mTranslationStem.at<double>(0, 0);
	translate[1] = mTranslationStem.at<double>(1, 0);
	translate[2] = mTranslationStem.at<double>(2, 0);

	itk::Rigid3DTransform<double>::Pointer transform = itk::VersorRigid3DTransform<double>::New();
	transform->SetMatrix(rotation);
	transform->SetOffset(translate);

	return transform;
}

cv::Mat HipFemurImplantsMatchInfo::Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer transform) const
{
	itk::Matrix< double, 3, 3 > rotation = transform->GetMatrix();
	double* matrix = new double[9];

	matrix[0] = rotation[0][0];
	matrix[1] = rotation[0][1];
	matrix[2] = rotation[0][2];

	matrix[3] = rotation[1][0];
	matrix[4] = rotation[1][1];
	matrix[5] = rotation[1][2];

	matrix[6] = rotation[2][0];
	matrix[7] = rotation[2][1];
	matrix[8] = rotation[2][2];

	cv::Mat matrixCV(3, 3, CV_64FC1, matrix);
	return matrixCV;
}

cv::Mat HipFemurImplantsMatchInfo::Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer transform) const
{
	itk::Vector< double, 3 > translate = transform->GetOffset();
	cv::Mat result(3, 1, CV_64FC1);

	result.at <double>(0, 0) = translate[0];
	result.at <double>(1, 0) = translate[1];
	result.at <double>(2, 0) = translate[2];
	return result;
}

Plane HipFemurImplantsMatchInfo::TransformPlane(const Plane& plane, const cv::Mat& rotation, const cv::Mat& translation) const
{
	cv::Mat transformNormalVector = rotation * plane.getNormalVectorMat();
	cv::Mat transformPoint = (rotation * plane.getPointMat()) + translation;
	Plane transformPlane;
	transformPlane.init(Point(transformNormalVector), Point(transformPoint));
	return transformPlane;
}

double HipFemurImplantsMatchInfo::getStemVersion() const
{
	cv::Mat neckAxisMat = mRotationStem * mImplantStem.getVectorNeckToHead().ToMatPoint();
	Point neckAxis = Point(neckAxisMat);

	return mPelvis.getFemurVersion(neckAxis);
}

double HipFemurImplantsMatchInfo::getCombinedOffsetDistance() const
{
	Plane coronal;
	Point temp = mImplantStem.getCanalAxisPoint() + 10. * mImplantStem.getVectorInfSup();
	coronal.init(mImplantStem.getHeadCenter(), mImplantStem.getCanalAxisPoint(), temp);
	coronal = TransformPlane(coronal, mRotationStem, mTranslationStem);

	Line femoralCanalAxes(mImplantStem.getVectorInfSup(), mImplantStem.getCanalAxisPoint());
	femoralCanalAxes.TransformLine(mRotationStem, mTranslationStem);

	return femoralCanalAxes.getDistanceFromPoint(coronal.getProjectionPoint(mPelvis.getPubicJoin()));
}

double HipFemurImplantsMatchInfo::getHipLengthDistance() const
{
	Plane coronal;
	Point temp = mImplantStem.getCanalAxisPoint() + 10. * mImplantStem.getVectorInfSup();
	coronal.init(mImplantStem.getHeadCenter(), mImplantStem.getCanalAxisPoint(), temp);
	coronal = TransformPlane(coronal, mRotationStem, mTranslationStem);

	Line lineASIS = Line::makeLineWithPoints(coronal.getProjectionPoint(mPelvis.getRightASIS()), coronal.getProjectionPoint(mPelvis.getLeftASIS()));
	Point ref = coronal.getProjectionPoint(mPelvis.getFemurOperationSide().getLesserTrochanter());

	return lineASIS.getDistanceFromPoint(ref);
}

