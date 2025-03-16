
#include "FemurImplantMatch.hpp"
#include <fstream>
#include "ImplantsException.hpp"
#include "ConvexHull.hpp"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkPlaneCollection.h"
#include "vtkClipClosedSurface.h"
#include "ImplantTools.hpp"

using namespace UKA::IMPLANTS;

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

void FemurImplantMatch::init(const FemurImplant& implantFemur, const Knee& knee)
{
	if (isInit == true)
	{
		throw ImplantExceptionCode::ALREADY_INITIALIZED_FEMUR_IMPLANT_MATCH;
	}
	this->implantFemur = implantFemur;
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

	Point implantFemurAxis = implantFemur.getDirectVectorFemurAxis();
	Point implantAP = implantFemur.getDirectVectorAP();
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
		tCenter = knee.getFemurKneeCenter() + 0.6 * (knee.getMedialEpicondylePerp() - knee.getFemurKneeCenter());
	}
	else
	{
		tCenter = knee.getFemurKneeCenter() + 0.5 * (knee.getLateralEpicondyle() - knee.getFemurKneeCenter());
	}

	/*
	if (knee.getSurgerySide() == SurgerySideEnum::KMedial)
	{
		tCenter = knee.getMedialEpicondyle();
	}
	else
	{
		tCenter = knee.getLateralEpicondyle();
	}
	*/

	//cv::Mat tempTuber = rotationMatrixTibiaImplant * implantTibia.getPointTuber().ToMatPoint() + translationMatrixTibiaImplant;
	//cv::Mat tempPCL = rotationMatrixTibiaImplant * implantTibia.getPointPCL().ToMatPoint() + translationMatrixTibiaImplant;
	//cv::Mat tempSide = rotationMatrixTibiaImplant * implantTibia.getPlateauRefPointDown().ToMatPoint() + translationMatrixTibiaImplant;

	//Plane baseTibiaImplant;
	//baseTibiaImplant.init(Point(tempTuber), Point(tempPCL), Point(tempSide));
	//Plane sidePlane = baseTibiaImplant.getPerpendicularPlane(Point(tempTuber), Point(tempPCL));
	//sidePlane.reverseByPoint(Point(tempSide));
	//sidePlane.movePlaneOnNormal((implantFemur.getWidthSize() / 2.) + 1.); //One is added to prevent the femur implant from adhering to the limit of the tibia implant.
	//tCenter = sidePlane.getProjectionPoint(Point(tempSide));

	cv::Mat kneeNormalVectorPlaneA = rotationMatrix * implantFemur.getPosterior().getNormalVectorMat();
	Plane kneePlaneA;
	kneePlaneA.init(Point(kneeNormalVectorPlaneA), knee.getMoveCondyle(implantFemur.getImplantInfo()));

	cv::Mat kneeNormalVectorMidPlane = rotationMatrix * implantFemur.getMidPlane().getNormalVectorMat();
	Plane kneeMidPlane;
	kneeMidPlane.init(Point(kneeNormalVectorMidPlane), tCenter);

	cv::Mat kneeNormalVectorPlaneC = rotationMatrix * implantFemur.getDistalPlane().getNormalVectorMat();
	Plane kneePlaneC;
	kneePlaneC.init(Point(kneeNormalVectorPlaneC), knee.getInferiorMoveFemurPoint(implantFemur.getImplantInfo()));

	cv::Mat pSeudoExpKneePointA = rotationMatrix * implantFemur.getPosterior().getPointMat();
	cv::Mat pSeudoExpKneeMidPoint = rotationMatrix * implantFemur.getMidPlane().getPointMat();
	cv::Mat pSeudoExpKneePointC = rotationMatrix * implantFemur.getDistalPlane().getPointMat();

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
	return result;

	/*
	cv::Mat normalRef = rotationMatrix * implantFemur.getMidPlane().getNormalVector().ToMatPoint();
	Plane planeSideRef;
	planeSideRef.init(Point(normalRef), tCenter);
	planeSideRef.reverseByPoint(knee.getFemurKneeCenter());

	double distanceToCenterKnee = planeSideRef.getDistanceFromPoint(knee.getFemurKneeCenter());
	double moveDist = distanceToCenterKnee - (implantFemur.getWidthSize());

	translationMatrix = translationMatrix + moveDist * planeSideRef.getNormalVector().ToMatPoint();

	return result;
	*/

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

Point FemurImplantMatch::getPointsOnPlane(const Plane& myPlane, std::vector<Point>& points) const
{
	Point cortex = knee.getCortexRef();
	Point hip = knee.getHipCenter();
	Point femurKnee = knee.getFemurKneeCenter();
	Point directVector = hip - femurKnee;
	Plane Transverse;
	Transverse.init(directVector, cortex);
	int transverseSign;
	if (Transverse.eval(hip) > 0)
		transverseSign = -1.0;
	else
		transverseSign = 1.0;

	Point midPoint, tempPoint;
	points.clear();

	vtkSmartPointer<vtkPolyData> contour = ImplantTools::getContours(knee.GetFemurPoly(), myPlane.getNormalVector(), myPlane.getPoint());
	vtkSmartPointer<vtkPoints> pointsCut = contour->GetPoints();
	int tSize = pointsCut->GetNumberOfPoints();
	for (int i = 0; i < tSize; i++)
	{
		double pnt[3];
		pointsCut->GetPoint(i, pnt);
		Point myPoint(pnt[0], pnt[1], pnt[2]);

		if (transverseSign * Transverse.eval(myPoint) > 0)
		{
			midPoint = midPoint + myPoint;
			points.push_back(myPoint);
		}
	}

	if (points.size() > 0)
	{
		midPoint = midPoint / double(points.size());
	}

	return midPoint;
}

Point FemurImplantMatch::getPointsOnPlane(const Plane& myPlane, std::vector<Point>& allPoints, std::vector<Point>& pointsLat, std::vector<Point>& pointsMed) const
{
	Point hip = knee.getHipCenter();
	Point femurKnee = knee.getFemurKneeCenter();
	Point cortex = femurKnee + 0.25 * (hip - femurKnee);
	Point directVector = hip - femurKnee;
	Plane Transverse, midPlane;
	Transverse.init(directVector, cortex);
	int transverseSign;
	if (Transverse.eval(hip) > 0)
		transverseSign = -1.0;
	else
		transverseSign = 1.0;

	midPlane.init(knee.getFemurVectorTEA(), femurKnee);
	midPlane.reverseByPoint(knee.getLateralEpicondyle());

	Point midPoint, tempPoint;
	pointsLat.clear();
	pointsMed.clear();

	auto allContours = ImplantTools::getAllContours(knee.GetFemurPoly(), myPlane.getNormalVector(), myPlane.getPoint());
	bool separated = true;

	for (int i = 0; i < allContours.size(); i++)
	{
		//ImplantTools::show(allContours[i].first, allContours[i].first);
		vtkSmartPointer<vtkPoints> pointsCut = allContours[i].second;
		int tSize = pointsCut->GetNumberOfPoints();
		int lateralCount = 0;
		int medialCount = 0;
		midPlane.countPositiveAndNegativePoints(pointsCut, lateralCount, medialCount);

		for (int j = 0; j < tSize; j++)
		{
			double pnt[3];
			pointsCut->GetPoint(j, pnt);
			Point myPoint(pnt[0], pnt[1], pnt[2]);

			if (transverseSign * Transverse.eval(myPoint) > 0)
			{
				midPoint = midPoint + myPoint;

				allPoints.push_back(myPoint);

				if (lateralCount > medialCount && separated == true)
				{
					float comp = medialCount * 2.5;
					if (lateralCount >= comp)
					{
						pointsLat.push_back(myPoint);
					}
					else
					{
						pointsLat.clear();
						pointsMed.clear();
						separated = false;
					}
				}
				else if (medialCount > lateralCount && separated == true)
				{
					float comp = lateralCount * 2.5;
					if (medialCount >= comp)
					{
						pointsMed.push_back(myPoint);
					}
					else
					{
						pointsLat.clear();
						pointsMed.clear();
						separated = false;
					}
				}

			}
		}

	}

	if (allPoints.size() > 0)
	{
		midPoint = midPoint / double(allPoints.size());
	}

	return midPoint;
}

std::vector<PointTypeITK> FemurImplantMatch::GetHullPoints(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut, BoneArea id, double distanceSide, double distanceTop, double angleLateral, double angleMedial, int amount) const
{
	std::vector<Point> projectedPoints;
	Point midPointPlane;
	std::vector<PointTypeITK> hull;
	double resizeVector = 1000000.0;

	Point centerP1, centerP2, topPoint, downPoint, lateralPoint, medialPoint, lateralSide, medialSide;
	Plane midPlane = finalTransformPlane(implantFemur.getMidPlane(), pTransformIn);
	Point fromPostToAntVector = knee.getFemurDirectVectorAP();

	lateralPoint = knee.getLateralEpicondyle();

	medialPoint = knee.getMedialEpicondylePerp();

	Point fromMedToLatVector = lateralPoint - medialPoint;
	fromMedToLatVector = fromMedToLatVector / sqrt(fromMedToLatVector.dot(fromMedToLatVector));

	Point anterior = knee.getFemurKneeCenter() + (resizeVector * fromPostToAntVector);

	Point posterior = knee.getFemurKneeCenter() - (resizeVector * fromPostToAntVector);

	lateralPoint = knee.getFemurKneeCenter() + (resizeVector * fromMedToLatVector);

	medialPoint = knee.getFemurKneeCenter() - (resizeVector * fromMedToLatVector);
	Plane currentPlane, sagitalAnatomicPlane;
	sagitalAnatomicPlane.init(knee.getFemurVectorTEA(), knee.getFemurKneeCenter());
	cv::Mat myRotation;
	Point myNormal, myNormalTemp;
	std::vector<Point> pointsLatTemp, pointsMedTemp;

	if (id == KPosteriorPlane)
	{
		currentPlane = finalTransformPlane(implantFemur.getPosterior(), pTransformIn);
		projectedPoints.clear();

		//midPointPlane = getPointsOnPlane(currentPlane, projectedPoints);
		midPointPlane = getPointsOnPlane(currentPlane, projectedPoints, pointsLatTemp, pointsMedTemp);

		myNormalTemp = knee.getFemurKneeCenter() - resizeVector * knee.getFemurDirectVectorAP();
		myNormal = currentPlane.getProjectionPoint(myNormalTemp) - myNormalTemp;

		currentPlane.fixNormalVector(myNormal);

		if (projectedPoints.size() > 15)
		{
			centerP1 = midPointPlane;
			centerP2 = currentPlane.getProjectionPoint(knee.getHipCenter());

			myRotation = getTransformToRobot(currentPlane, sagitalAnatomicPlane, centerP1, centerP2);

			centerP1 = midPlane.getProjectionPoint(centerP1);
			centerP2 = midPlane.getProjectionPoint(centerP2);

			topPoint = currentPlane.getProjectionPoint(knee.getHipCenter());
			downPoint = currentPlane.getProjectionPoint(knee.getAnkleCenter());
			lateralSide = currentPlane.getProjectionPoint(lateralPoint);
			medialSide = currentPlane.getProjectionPoint(medialPoint);

		}
		else
		{
			return hull;
		}
	}
	else
	{
		currentPlane = finalTransformPlane(implantFemur.getBestPlaneToCurvePoints(), pTransformIn);
		projectedPoints.clear();
		midPointPlane = getPointsOnPlane(currentPlane, projectedPoints);

		myNormalTemp = knee.getAnkleCenter();
		myNormal = myNormalTemp - currentPlane.getProjectionPoint(myNormalTemp);

		currentPlane.fixNormalVector(myNormal);

		if (projectedPoints.size() > 10)
		{
			centerP1 = midPointPlane;
			centerP2 = knee.getFemurKneeCenter() - resizeVector * knee.getFemurDirectVectorAP();
			centerP2 = currentPlane.getProjectionPoint(centerP2);

			myRotation = getTransformToRobot(currentPlane, sagitalAnatomicPlane, centerP1, centerP2);
		}
		else
		{
			return hull;
		}
	}

	std::vector<Point> vertices, cutPoints;

	if (angleLateral > 45)
	{
		angleLateral = 45;
	}

	if (angleLateral < 0)
	{
		angleLateral = 0;
	}

	if (angleMedial > 45)
	{
		angleMedial = 45;
	}

	if (angleMedial < 0)
	{
		angleMedial = 0;
	}

	double angleLatRad = ((90.0 - angleLateral) * PI) / 180.0;
	double angleMedRad = ((90.0 - angleMedial) * PI) / 180.0;

	if (id == KPosteriorPlane)
	{
		if (knee.getSurgerySide() == SurgerySideEnum::KLateral && pointsLatTemp.size() > 3)
		{
			getCurveLikeU(pointsLatTemp, downPoint, lateralSide, medialSide, topPoint, midPlane, currentPlane, myRotation, vertices, distanceSide, distanceTop, amount);
		}
		else if (knee.getSurgerySide() == SurgerySideEnum::KMedial && pointsMedTemp.size() > 3)
		{
			getCurveLikeU(pointsMedTemp, downPoint, lateralSide, medialSide, topPoint, midPlane, currentPlane, myRotation, vertices, distanceSide, distanceTop, amount);
		}
		else
		{
			getCurveLikeU(pointsLatTemp, downPoint, lateralSide, medialSide, topPoint, midPlane, currentPlane, myRotation, vertices, distanceSide, distanceTop, amount);
		}
	}
	else
	{
		cv::Mat rotation = Rigid3DTransformToCVRotation(pTransformIn);
		cv::Mat translation = Rigid3DTransformToCVTranslation(pTransformIn);
		// No esta actualizado aun en la nube de china
		vertices = ConvexHull::interpolateSpline(implantFemur.getAllSidePointsInOrder(rotation, translation, distanceSide), amount);
	}

	Point initExtreme, endExtreme;

	if (vertices.size() > 1)
	{
		initExtreme = vertices[0];
		endExtreme = vertices[vertices.size() - 1];

		if (id == KAnteriorAndDistalCurve)
		{
			initExtreme = currentPlane.getProjectionPoint(initExtreme);
			endExtreme = currentPlane.getProjectionPoint(endExtreme);
		}

		cv::Mat initExtremeMat, endExtremeMat;

		initExtremeMat = myRotation * initExtreme.ToMatPoint();
		endExtremeMat = myRotation * endExtreme.ToMatPoint();

		initExtreme = Point(initExtremeMat);
		endExtreme = Point(endExtremeMat);

		if (endExtreme.y < initExtreme.y)
		{
			std::reverse(vertices.begin(), vertices.end());
		}
	}

	hull = increaseVectorToAmount(vertices, amount);

	itk::Vector< double, 3 > translate;

	//ImplantTools::fitEllipse(vertices, currentPlane.getNormalVector(), midPointPlane);

	midPointPlane = ImplantTools::getPolygonCenter(vertices, currentPlane.getNormalVector());

	/*auto vectorTest = vertices;
	vectorTest.push_back(midPointPlane);
	ImplantTools::show(knee.GetFemurPoly(), vectorTest);*/

	if (projectedPoints.size() > 0)
	{
		midPointPlane = currentPlane.getProjectionPoint(midPointPlane);
		Point tTemp = myRotation * (midPointPlane.ToMatPoint());
		translate[0] = -tTemp.x;
		translate[1] = -tTemp.y;
		translate[2] = -tTemp.z;
	}
	else
	{
		translate[0] = 0;
		translate[1] = 0;
		translate[2] = 0;
	}

	itk::Matrix< double, 3, 3 > rotationITK;
	rotationITK[0][0] = myRotation.at<double>(0, 0);
	rotationITK[0][1] = myRotation.at<double>(0, 1);
	rotationITK[0][2] = myRotation.at<double>(0, 2);

	rotationITK[1][0] = myRotation.at<double>(1, 0);
	rotationITK[1][1] = myRotation.at<double>(1, 1);
	rotationITK[1][2] = myRotation.at<double>(1, 2);

	rotationITK[2][0] = myRotation.at<double>(2, 0);
	rotationITK[2][1] = myRotation.at<double>(2, 1);
	rotationITK[2][2] = myRotation.at<double>(2, 2);

	pTransformOut->SetMatrix(rotationITK);
	pTransformOut->SetOffset(translate);

	return hull;
}

Plane FemurImplantMatch::finalTransformPlane(const Plane& plane, const itk::Rigid3DTransform<>::Pointer pTransform) const
{
	cv::Mat rotation = Rigid3DTransformToCVRotation(pTransform);
	cv::Mat translation = Rigid3DTransformToCVTranslation(pTransform);

	cv::Mat transformNormalVector = rotation * plane.getNormalVectorMat();
	cv::Mat transformPoint = (rotation * plane.getPointMat()) + translation;
	Plane transformPlane;
	transformPlane.init(Point(transformNormalVector), Point(transformPoint));
	return transformPlane;
}

void FemurImplantMatch::getCurveLikeU(const std::vector<Point>& points, const Point& downPoint, const Point& lateralPoint, const Point& medialPoint, const Point& topPoint, const Plane& midPlane, const Plane& currentPlane, const cv::Mat& pRotation, std::vector<Point>& vertices, double distanceSide, double distanceTop, int amount) const
{
	ConvexHullFeatures hullFeatures = getIncreaseBorder(points, downPoint, lateralPoint, medialPoint, topPoint, midPlane, currentPlane, pRotation, distanceSide, distanceSide, distanceTop);

	vertices = ConvexHull::interpolateSpline(hullFeatures.convexHull, amount);
}

FemurImplantMatch::ConvexHullFeatures FemurImplantMatch::getIncreaseBorder(const std::vector<Point>& points, const Point& downPoint, const Point& lateralPoint, const Point& medialPoint, const Point& topPoint, const Plane& midPlane, const Plane& currentPlane, const cv::Mat& pRotation, double distanceSideLat, double distanceSideMed, double distanceTop, double downLatCornerOut, double downMedCornerOut) const
{
	double maxDist = 15.;

	if (abs(distanceSideLat) > maxDist)
	{
		distanceSideLat = maxDist;
	}

	if (abs(distanceSideMed) > maxDist)
	{
		distanceSideMed = maxDist;
	}

	if (abs(distanceTop) > maxDist)
	{
		distanceTop = maxDist;
	}

	Point vectorAP = midPlane.getProjectionPoint(topPoint) - midPlane.getProjectionPoint(downPoint);
	vectorAP.normalice();

	Line topLine(midPlane.getNormalVector(), topPoint);
	Line downLine(midPlane.getNormalVector(), downPoint);
	Line medialLine(vectorAP, medialPoint);
	Line lateralLine(vectorAP, lateralPoint);

	double topDistance = -1.0;
	double downDistance = -1.0;
	double medialDistance = -1.0;
	double lateralDistance = -1.0;

	double topTemp, downTemp, lateralTemp, medialTemp;
	Point topPointConv, downPointConv, latPointConv, medPointConv, centerPoint;

	ConvexHull allHull(points, pRotation);
	std::vector<Point> fullConvex = allHull.GetConvexHull();
	int tSize = fullConvex.size();

	if (tSize < 3)
	{
		throw ImplantExceptionCode::CAN_NOT_DETERMINE_CONVEX_HULL;
	}

	std::vector<Point>::const_iterator it1, it2;
	it1 = fullConvex.begin();
	it2 = fullConvex.end();

	for (; it1 != it2; ++it1)
	{
		topTemp = topLine.getDistanceFromPoint(*it1);
		downTemp = downLine.getDistanceFromPoint(*it1);
		lateralTemp = lateralLine.getDistanceFromPoint(*it1);
		medialTemp = medialLine.getDistanceFromPoint(*it1);

		if (topTemp < topDistance || topDistance < 0)
		{
			topDistance = topTemp;
			topPointConv = *it1;
		}

		if (downTemp < downDistance || downDistance < 0)
		{
			downDistance = downTemp;
			downPointConv = *it1;
		}

		if (lateralTemp < lateralDistance || lateralDistance < 0)
		{
			lateralDistance = lateralTemp;
			latPointConv = *it1;
		}

		if (medialTemp < medialDistance || medialDistance < 0)
		{
			medialDistance = medialTemp;
			medPointConv = *it1;
		}

		centerPoint = centerPoint + (*it1);
	}

	topLine.setPoint(topPointConv);
	downLine.setPoint(downPointConv);
	lateralLine.setPoint(latPointConv);
	medialLine.setPoint(medPointConv);
	centerPoint = centerPoint / double(tSize);

	int topLatCorner = ImplantTools::GetCornerPointOnContour(fullConvex, centerPoint, topLine.getProjectPoint(centerPoint), lateralLine.getProjectPoint(centerPoint));
	int topMedCorner = ImplantTools::GetCornerPointOnContour(fullConvex, centerPoint, topLine.getProjectPoint(centerPoint), medialLine.getProjectPoint(centerPoint));
	int downLatCorner = ImplantTools::GetCornerPointOnContour(fullConvex, centerPoint, downLine.getProjectPoint(centerPoint), lateralLine.getProjectPoint(centerPoint));
	int downMedCorner = ImplantTools::GetCornerPointOnContour(fullConvex, centerPoint, downLine.getProjectPoint(centerPoint), medialLine.getProjectPoint(centerPoint));

	////////////////////////////////////////////////////////////////////////////////////
	//auto poly1 = ImplantTools::getContours(knee.GetFemurPoly(), currentPlane.getNormalVector(), currentPlane.getPoint());
	//ImplantTools::show(poly1, fullConvex, true);
	////////////////////////////////////////////////////////////////////////////////////

	if (topLatCorner < 0 || topMedCorner < 0 || downLatCorner < 0 || downMedCorner < 0)
	{
		throw ImplantExceptionCode::CAN_NOT_DETERMINE_CONVEX_HULL_CORNERS;
	}

	ConvexHullFeatures result;
	result.centerPoint = centerPoint;
	result.lateralTopPoint = fullConvex[topLatCorner];
	result.medialTopPoint = fullConvex[topMedCorner];

	Plane obliqueLatTopUp, obliqueMedTopUp, obliqueLatTopSide, obliqueLatDownSide, obliqueMedTopSide, obliqueMedDownSide;

	obliqueLatTopUp = currentPlane.getPerpendicularPlane(centerPoint, fullConvex[topLatCorner]);
	obliqueMedTopUp = currentPlane.getPerpendicularPlane(centerPoint, fullConvex[topMedCorner]);
	obliqueLatTopUp.reverseByPoint(fullConvex[topMedCorner]);
	obliqueMedTopUp.reverseByPoint(fullConvex[topLatCorner]);

	obliqueLatTopSide = currentPlane.getPerpendicularPlane(centerPoint, fullConvex[topLatCorner]);
	obliqueLatDownSide = currentPlane.getPerpendicularPlane(centerPoint, fullConvex[downLatCorner]);
	obliqueLatTopSide.reverseByPoint(fullConvex[downLatCorner]);
	obliqueLatDownSide.reverseByPoint(fullConvex[topLatCorner]);

	obliqueMedTopSide = currentPlane.getPerpendicularPlane(centerPoint, fullConvex[topMedCorner]);
	obliqueMedDownSide = currentPlane.getPerpendicularPlane(centerPoint, fullConvex[downMedCorner]);
	obliqueMedTopSide.reverseByPoint(fullConvex[downMedCorner]);
	obliqueMedDownSide.reverseByPoint(fullConvex[topMedCorner]);

	Plane myMidPlane = midPlane;
	myMidPlane.movePlane(medPointConv);
	myMidPlane.reverseByPoint(latPointConv);
	myMidPlane.movePlane(centerPoint);

	for (int i = 0; i < tSize; i++)
	{
		if (i == topLatCorner || i == topMedCorner)
		{
			fullConvex[i] = fullConvex[i] + distanceTop * vectorAP;

			if (i == topLatCorner)
			{
				fullConvex[i] = fullConvex[i] + distanceSideLat * myMidPlane.getNormalVector();
			}
			else
			{
				fullConvex[i] = fullConvex[i] - distanceSideMed * myMidPlane.getNormalVector();
			}
			continue;
		}

		if (i == downLatCorner || i == downMedCorner)
		{
			if (i == downLatCorner)
			{
				fullConvex[i] = fullConvex[i] + distanceSideLat * myMidPlane.getNormalVector();
				if (downLatCornerOut >= 0)
				{
					Point temp = downLine.getProjectPoint(fullConvex[i]);
					temp = lateralLine.getProjectPoint(temp);
					temp = temp + downLatCornerOut * distanceSideLat * myMidPlane.getNormalVector();
					fullConvex[i] = temp;
				}
			}
			else
			{
				fullConvex[i] = fullConvex[i] - distanceSideMed * myMidPlane.getNormalVector();
				if (downMedCornerOut >= 0)
				{
					Point temp = downLine.getProjectPoint(fullConvex[i]);
					temp = medialLine.getProjectPoint(temp);
					temp = temp - downMedCornerOut * distanceSideMed * myMidPlane.getNormalVector();
					fullConvex[i] = temp;
				}
			}
			continue;
		}


		if (obliqueLatTopUp.eval(fullConvex[i]) > 0 && obliqueMedTopUp.eval(fullConvex[i]) > 0)
		{
			fullConvex[i] = fullConvex[i] + distanceTop * vectorAP;
		}
		else if (obliqueLatTopSide.eval(fullConvex[i]) > 0 && obliqueLatDownSide.eval(fullConvex[i]) > 0)
		{
			fullConvex[i] = fullConvex[i] + distanceSideLat * myMidPlane.getNormalVector();
		}
		else if (obliqueMedTopSide.eval(fullConvex[i]) > 0 && obliqueMedDownSide.eval(fullConvex[i]) > 0)
		{
			fullConvex[i] = fullConvex[i] - distanceSideMed * myMidPlane.getNormalVector();
		}
	}

	ConvexHull convHull(fullConvex, pRotation);
	fullConvex = convHull.GetConvexHull();
	tSize = fullConvex.size();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//auto poly2 = ImplantTools::getContours(knee.GetFemurPoly(), currentPlane.getNormalVector(), currentPlane.getPoint());
	//ImplantTools::show(poly2, convHull.getChangeDirectionPoints());
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<Point> changeDirPoints = convHull.getChangeDirectionPoints();

	topLine.setPoint((topLine.getPoint() + distanceTop * vectorAP));
	lateralLine.setPoint((lateralLine.getPoint() + distanceSideLat * myMidPlane.getNormalVector()));
	medialLine.setPoint((medialLine.getPoint() - distanceSideMed * myMidPlane.getNormalVector()));

	if (tSize == 0)
	{
		throw ImplantExceptionCode::CAN_NOT_DETERMINE_CONVEX_HULL;
	}

	downLatCorner = ImplantTools::GetCornerPointOnContour(fullConvex, centerPoint, downLine.getProjectPoint(centerPoint), lateralLine.getProjectPoint(centerPoint));
	downMedCorner = ImplantTools::GetCornerPointOnContour(fullConvex, centerPoint, downLine.getProjectPoint(centerPoint), medialLine.getProjectPoint(centerPoint));

	if (downLatCorner < 0 || downMedCorner < 0)
	{
		throw ImplantExceptionCode::CAN_NOT_DETERMINE_CONVEX_HULL;
	}

	Plane obliqueLat = currentPlane.getPerpendicularPlane(centerPoint, fullConvex[downLatCorner]);
	Plane obliqueMed = currentPlane.getPerpendicularPlane(centerPoint, fullConvex[downMedCorner]);

	obliqueLat.reverseByPoint(fullConvex[downMedCorner]);
	obliqueMed.reverseByPoint(fullConvex[downLatCorner]);

	std::vector<Point> finalHull;
	int posLat, posMed;

	for (int i = 0; i < tSize; i++)
	{
		if (i == downLatCorner)
		{
			posLat = finalHull.size();
			result.lateralDownPos = posLat;
			finalHull.push_back(fullConvex[i]);
			continue;
		}
		else if (i == downMedCorner)
		{
			posMed = finalHull.size();
			result.medialDownPos = posMed;
			finalHull.push_back(fullConvex[i]);
			continue;
		}

		if (!(obliqueLat.eval(fullConvex[i]) > 0 && obliqueMed.eval(fullConvex[i]) > 0))
		{
			finalHull.push_back(fullConvex[i]);
		}
	}

	//auto poly = ImplantTools::getContours(knee.GetFemurPoly(), currentPlane.getNormalVector(), currentPlane.getPoint());
	//std::vector<Point> tempTest = { result.medialTopPoint, finalHull[result.medialDownPos] };
	//ImplantTools::show(poly, tempTest);

	if (abs(posMed - posLat) == 1)
	{
		if (posMed > posLat)
		{
			std::rotate(finalHull.begin(), finalHull.begin() + posMed, finalHull.end());

			result.lateralDownPos = ((result.lateralDownPos - posMed) >= 0) ? result.lateralDownPos - posMed : finalHull.size() + result.lateralDownPos - posMed;
			result.medialDownPos = ((result.medialDownPos - posMed) >= 0) ? result.medialDownPos - posMed : finalHull.size() + result.medialDownPos - posMed;
		}
		else
		{
			std::rotate(finalHull.begin(), finalHull.begin() + posLat, finalHull.end());

			result.lateralDownPos = ((result.lateralDownPos - posLat) >= 0) ? result.lateralDownPos - posLat : finalHull.size() + result.lateralDownPos - posLat;
			result.medialDownPos = ((result.medialDownPos - posLat) >= 0) ? result.medialDownPos - posLat : finalHull.size() + result.medialDownPos - posLat;
		}
	}

	//auto poly1 = ImplantTools::getContours(knee.GetFemurPoly(), currentPlane.getNormalVector(), currentPlane.getPoint());
	//std::vector<Point> tempTest1 = { result.medialTopPoint, finalHull[result.medialDownPos] };
	//ImplantTools::show(poly1, tempTest1);

	result.convexHull = finalHull;
	result.downLine = new Line(downLine.getDirectVector(), downLine.getPoint());
	result.topLine = new Line(topLine.getDirectVector(), topLine.getPoint());
	result.lateralLine = new Line(lateralLine.getDirectVector(), lateralLine.getPoint());
	result.medialLine = new Line(medialLine.getDirectVector(), medialLine.getPoint());

	return result;
}

std::vector<PointTypeITK> FemurImplantMatch::increaseVectorToAmount(const std::vector<Point>& points, int amount) const
{
	std::vector<PointTypeITK> result;
	if (points.size() >= amount || points.size() <= 1 || amount < 3)
	{
		auto it1 = points.begin();
		auto it2 = points.end();

		for (; it1 != it2; ++it1)
		{
			result.push_back((*it1).ToITKPoint());
		}
		return result;
	}

	int intervals = points.size() - 2;
	if (intervals == 0)
	{
		intervals = 1;
	}
	int stillPoints = amount - points.size();
	int interAmount = stillPoints / intervals;
	int lastPos = 0;
	int rest = stillPoints % intervals;
	Point a, b, c;
	double coef = 1.0 / double(interAmount + 1);
	double coefRest = 1.0 / double(rest + 1);
	for (int i = 0; i < points.size() - 1; i++)
	{
		lastPos = i;
		a = points[i];
		b = points[i + 1];
		result.push_back(a.ToITKPoint());
		for (int j = 1; j <= interAmount; j++)
		{
			c = a + double(j) * coef * (b - a);
			result.push_back(c.ToITKPoint());
		}
		intervals--;
		if (intervals == 0)
		{
			break;
		}
	}

	for (int i = lastPos + 1; i < points.size() - 1; i++)
	{
		if (i == lastPos + 1)
		{
			a = points[i];
			b = points[i + 1];

			result.push_back(a.ToITKPoint());
			for (int j = 1; j <= rest; j++)
			{
				c = a + double(j) * coefRest * (b - a);
				result.push_back(c.ToITKPoint());
			}
		}
		else
		{
			result.push_back((points[i]).ToITKPoint());
		}
	}
	result.push_back((points[points.size() - 1]).ToITKPoint());
	return result;
}

Point FemurImplantMatch::movePointAtNormal(const Point& movePoint, const Point& nextPoint, const cv::Mat& rotationZ, double distance, bool changeMove, bool clockWise) const
{
	cv::Mat moveMat = rotationZ * movePoint.ToMatPoint();
	cv::Mat nextMat = rotationZ * nextPoint.ToMatPoint();

	Point moveP = Point(moveMat);
	Point nextP = Point(nextMat);
	Point vector;
	if (clockWise == true)
	{
		vector = nextP - moveP;
	}
	else
	{
		vector = moveP - nextP;
	}

	cv::Point2d perpendicular2d(vector.y, -vector.x);
	perpendicular2d = perpendicular2d / sqrt(perpendicular2d.dot(perpendicular2d));
	Point perpendicular(perpendicular2d.x, perpendicular2d.y, 0);
	Point finalMove;
	if (changeMove == true)
	{
		finalMove = nextP + distance * perpendicular;
	}
	else
	{
		finalMove = moveP + distance * perpendicular;
	}

	cv::Mat resultMat = rotationZ.inv() * finalMove.ToMatPoint();
	return Point(resultMat);
}

