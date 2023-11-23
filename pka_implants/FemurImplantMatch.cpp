
#include "FemurImplantMatch.hpp"
#include <fstream>
#include "ImplantsException.hpp"
#include "ConvexHull.hpp"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkPlaneCollection.h"
#include "vtkClipClosedSurface.h"
#include "ImplantTools.hpp"

using namespace PKA::IMPLANTS;

inline cv::Mat Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer pTransform)
{
	itk::Matrix< double, 3, 3 > rotation = pTransform->GetMatrix();
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

inline cv::Mat Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer pTransform)
{
	itk::Vector< double, 3 > translate = pTransform->GetOffset();
	cv::Mat result(3, 1, CV_64FC1);

	result.at <double>(0, 0) = translate[0];
	result.at <double>(1, 0) = translate[1];
	result.at <double>(2, 0) = translate[2];
	return result;
}

inline cv::Mat GetRotateTransformAxisZ(Plane myPlane)
{
	//myPlane.normalizeNormalVector();
	Point normalXY(0, 0, 1);
	Point rotationAxis = normalXY.cross(myPlane.getNormalVector());
	rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
	double rotationAngle = Line::getAngleBetweenVectors(normalXY, myPlane.getNormalVector());
	cv::Mat rotation_1 = Line::getRotateMatrix(rotationAxis, -rotationAngle);
	cv::Mat rotation_2 = Line::getRotateMatrix(rotationAxis, rotationAngle);
	cv::Mat rotateVector_1 = rotation_1 * myPlane.getNormalVectorMat();
	cv::Mat rotateVector_2 = rotation_2 * myPlane.getNormalVectorMat();

	cv::Mat rotate;

	double distance_1 = Line::getAngleBetweenVectors(Point(rotateVector_1), normalXY);
	double distance_2 = Line::getAngleBetweenVectors(Point(rotateVector_2), normalXY);

	if (distance_1 < distance_2)
	{
		rotate = rotation_1;
	}
	else
	{
		rotate = rotation_2;
	}

	return rotate;
}

inline cv::Mat GetRotateTransformAxisX(Point sagitalVector, const cv::Mat& rotationZ)
{
	sagitalVector.normalice();
	cv::Mat finalSagitalMat = rotationZ * sagitalVector.ToMatPoint();
	Point finalSagitalVector = Point(finalSagitalMat);
	Point normalXY(1, 0, 0);
	Point rotationAxis = normalXY.cross(finalSagitalVector);
	rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
	double rotationAngle = Line::getAngleBetweenVectors(normalXY, finalSagitalVector);
	cv::Mat rotation_1 = Line::getRotateMatrix(rotationAxis, -rotationAngle);
	cv::Mat rotation_2 = Line::getRotateMatrix(rotationAxis, rotationAngle);
	cv::Mat rotateVector_1 = rotation_1 * finalSagitalMat;
	cv::Mat rotateVector_2 = rotation_2 * finalSagitalMat;

	cv::Mat rotate;

	double distance_1 = Line::getAngleBetweenVectors(Point(rotateVector_1), normalXY);
	double distance_2 = Line::getAngleBetweenVectors(Point(rotateVector_2), normalXY);

	if (distance_1 < distance_2)
	{
		rotate = rotation_1;
	}
	else
	{
		rotate = rotation_2;
	}

	return rotate;
}
inline cv::Mat getTransformToRobot(Plane currentPlane, const Plane& sagitalAnatomicPlane, const Point& P1, const Point& P2)
{
	Point sagitalP1 = sagitalAnatomicPlane.getProjectionPoint(P1);
	Point sagitalP2 = sagitalAnatomicPlane.getProjectionPoint(P2);
	sagitalP1 = currentPlane.getProjectionPoint(sagitalP1);
	sagitalP2 = currentPlane.getProjectionPoint(sagitalP2);
	Point sagitalVector = sagitalP2 - sagitalP1;

	//Se rota el plano hacia el eje Z.
	cv::Mat rotateZ = GetRotateTransformAxisZ(currentPlane);

	//Una vez el plano es perpendicular al eje Z, se gira hasta que su vector sagital sea paralelo al eje X.

	cv::Mat rotateX = GetRotateTransformAxisX(sagitalVector, rotateZ);
	cv::Mat finalMatrix = rotateX * rotateZ;
	return finalMatrix;
}

FemurImplantMatch::FemurImplantMatch()
{
	isInit = false;
}
void FemurImplantMatch::init(const FemurImplant& implant, const Knee& knee)
{
	if (isInit == true)
	{
		throw ImplantExceptionCode::ALREADY_INITIALIZED_FEMUR_IMPLANT_MATCH;
	}
	this->implant = implant;
	this->knee = knee;
	getRotationMatrix();
	bool result = getTranslationMatrix();
	if (result == false)
	{
		throw ImplantExceptionCode::FAILED_TRANSFORMATION_MATCH;
	}
	isInit = true;
}

void FemurImplantMatch::getRotationMatrix()
{
	std::vector<cv::Point3d> implantVectors;
	std::vector<cv::Point3d> kneeVectors;

	Point implantFemurAxis = implant.getDirectVectorFemurAxis();
	Point implantAP = implant.getDirectVectorAP();
	Point implantTEA = implantFemurAxis.cross(implantAP);
	implantTEA.normalice();

	implantVectors.push_back(implantTEA.ToCVPoint());
	implantVectors.push_back(implantFemurAxis.ToCVPoint());
	implantVectors.push_back(implantAP.ToCVPoint());

	Point kneeFemurAxis = knee.getDirectVectorFemurAxis();
	Point kneeAP = knee.getFemurDirectVectorAP();
	Point kneeTEA = kneeFemurAxis.cross(kneeAP);
	kneeTEA.normalice();

	kneeVectors.push_back(kneeTEA.ToCVPoint());
	kneeVectors.push_back(kneeFemurAxis.ToCVPoint());
	kneeVectors.push_back(kneeAP.ToCVPoint());

	cv::Mat implantMatrix = cv::Mat(implantVectors.size(), 3, CV_64F, implantVectors.data());
	cv::Mat kneeMatrix = cv::Mat(kneeVectors.size(), 3, CV_64F, kneeVectors.data());

	cv::Mat inverse = (implantMatrix.t()).inv();
	rotationMatrix = (kneeMatrix.t()) * inverse;
	//std::cout << "Rotation: " << rotationMatrix << std::endl;
}

bool FemurImplantMatch::getTranslationMatrix()
{
	Point tCenter;

	if (knee.getSurgerySide() == SurgerySideEnum::KMedial)
	{
		tCenter = knee.getMedialEpicondyle();
	}
	else
	{
		tCenter = knee.getLateralEpicondyle();
	}

	cv::Mat kneeNormalVectorPlaneA = rotationMatrix * implant.getPosterior().getNormalVectorMat();
	Plane kneePlaneA;
	kneePlaneA.init(Point(kneeNormalVectorPlaneA), knee.getMoveCondyle(implant.getImplantInfo()));

	cv::Mat kneeNormalVectorMidPlane = rotationMatrix * implant.getMidPlane().getNormalVectorMat();
	Plane kneeMidPlane;
	kneeMidPlane.init(Point(kneeNormalVectorMidPlane), tCenter);

	cv::Mat kneeNormalVectorPlaneC = rotationMatrix * implant.getDistalPlane().getNormalVectorMat();
	Plane kneePlaneC;
	kneePlaneC.init(Point(kneeNormalVectorPlaneC), knee.getInferiorMoveFemurPoint(implant.getImplantInfo()));

	cv::Mat pSeudoExpKneePointA = rotationMatrix * implant.getPosterior().getPointMat();
	cv::Mat pSeudoExpKneeMidPoint = rotationMatrix * implant.getMidPlane().getPointMat();
	cv::Mat pSeudoExpKneePointC = rotationMatrix * implant.getDistalPlane().getPointMat();

	Point pSeudoKneePointA(pSeudoExpKneePointA);
	Point pSeudoKneeMidPoint(pSeudoExpKneeMidPoint);
	Point pSeudoKneePointC(pSeudoExpKneePointC);

	std::vector<cv::Point3d> normalVectors;
	std::vector<double> biasVector;
	double bias = 0.0;
	normalVectors.push_back(kneePlaneA.getNormalVector().ToCVPoint());
	normalVectors.push_back(kneeMidPlane.getNormalVector().ToCVPoint());
	normalVectors.push_back(kneePlaneC.getNormalVector().ToCVPoint());
	cv::Mat A(normalVectors.size(), 3, CV_64F, normalVectors.data());
	bias = -(kneePlaneA.getBias() + (pSeudoKneePointA.dot(kneePlaneA.getNormalVector())));
	biasVector.push_back(bias);
	bias = -(kneeMidPlane.getBias() + (pSeudoKneeMidPoint.dot(kneeMidPlane.getNormalVector())));
	biasVector.push_back(bias);
	bias = -(kneePlaneC.getBias() + (pSeudoKneePointC.dot(kneePlaneC.getNormalVector())));
	biasVector.push_back(bias);
	cv::Mat B(biasVector.size(), 1, CV_64F, biasVector.data());
	bool result = cv::solve(A, B, translationMatrix);

	cv::Mat normalRef = rotationMatrix * implant.getMidPlane().getNormalVector().ToMatPoint();
	Plane planeSideRef;
	planeSideRef.init(Point(normalRef), tCenter);
	planeSideRef.reverseByPoint(knee.getFemurKneeCenter());

	double distanceToCenterKnee = planeSideRef.getDistanceFromPoint(knee.getFemurKneeCenter());
	//double moveDist = distanceToCenterKnee - (implant.getWidthSize() / 2.);
	double moveDist = distanceToCenterKnee - (implant.getWidthSize()); // usando esta distancia temporal ****************

	translationMatrix = translationMatrix + moveDist * planeSideRef.getNormalVector().ToMatPoint();

	return result;

	/////////////////////////////////////////////////////////////// Temporal para probar *******************

	/*if (knee.getSurgerySide() == SurgerySideEnum::KMedial)
	{
		tCenter = knee.getMedialInferiorFemurPoint();
	}
	else
	{
		tCenter = knee.getLateralInferiorFemurPoint();
	}

	translationMatrix = tCenter.ToMatPoint() - rotationMatrix * implant.getRodBasePoint().ToMatPoint();

	return true;*/
}

itk::Matrix< double, 3, 3 > FemurImplantMatch::GetRotationMatrix() const
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

itk::Vector< double, 3 > FemurImplantMatch::GetTranslationMatrix() const
{
	itk::Vector< double, 3 > translation;
	translation[0] = translationMatrix.at <double>(0, 0);
	translation[1] = translationMatrix.at <double>(1, 0);
	translation[2] = translationMatrix.at <double>(2, 0);

	return translation;
}
