#include "TibiaSpacerImplantMatch.hpp"
#include "ImplantsException.hpp"
#include "ImplantTools.hpp"

using namespace UKA::IMPLANTS;

TibiaSpacerImplantMatch::TibiaSpacerImplantMatch()
{
	isInit = false;
}

void TibiaSpacerImplantMatch::init(const TibiaImplant& tibiaImplant, const TibiaSpacerImplant& spacerImplant)
{
	if (isInit == true)
	{
		throw ImplantExceptionCode::ALREADY_INITIALIZED_TIBIA_SPACER_IMPLANT_MATCH;
	}
	this->implant = tibiaImplant;
	this->spacerImplant = spacerImplant;

	makeRotationMatrix();
	makeTranslationMatrix();
	isInit = true;
}

void TibiaSpacerImplantMatch::makeRotationMatrix()
{
	std::vector<cv::Point3d> spacerVectors;
	std::vector<cv::Point3d> tibiaVectors;

	Point implantSpacerNormal = spacerImplant.getSpacerNormalVector();
	Point implantSpacerAP = spacerImplant.getSpacerVectorAP();
	Point implantSpacerCross = implantSpacerNormal.cross(implantSpacerAP);
	spacerVectors.push_back(implantSpacerNormal.ToCVPoint());
	spacerVectors.push_back(implantSpacerAP.ToCVPoint());
	spacerVectors.push_back(implantSpacerCross.ToCVPoint());

	Point implantTibiaNormal = implant.getTibiaNormalVector();
	Point implantTibiaAP = implant.getTibiaVectorAP();
	Point implantTibiaCross = implantTibiaNormal.cross(implantTibiaAP);

	tibiaVectors.push_back(implantTibiaNormal.ToCVPoint());
	tibiaVectors.push_back(implantTibiaAP.ToCVPoint());
	tibiaVectors.push_back(implantTibiaCross.ToCVPoint());

	cv::Mat implantSpacerMatrix = cv::Mat(spacerVectors.size(), 3, CV_64F, spacerVectors.data());
	cv::Mat implantTibiaMatrix = cv::Mat(tibiaVectors.size(), 3, CV_64F, tibiaVectors.data());
	cv::Mat inverse = (implantSpacerMatrix.t()).inv();
	rotationMatrix = (implantTibiaMatrix.t()) * inverse;
}

void TibiaSpacerImplantMatch::makeTranslationMatrix()
{
	translationMatrix = implant.getCentralPointUp().ToMatPoint() - (rotationMatrix * spacerImplant.getSpacerKneeCenter().ToMatPoint());
}

itk::Matrix< double, 3, 3 > TibiaSpacerImplantMatch::GetRotationMatrix() const
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

itk::Vector< double, 3 > TibiaSpacerImplantMatch::GetTranslationMatrix() const
{
	itk::Vector< double, 3 > translation;
	translation[0] = translationMatrix.at <double>(0, 0);
	translation[1] = translationMatrix.at <double>(1, 0);
	translation[2] = translationMatrix.at <double>(2, 0);

	//Point tt = Point(translationMatrix);
	//tt = tt + 3.0 * knee.getTibiaVectorLateralTEA();

	//translation[0] = tt.x;
	//translation[1] = tt.y;
	//translation[2] = tt.z;

	return translation;
}

