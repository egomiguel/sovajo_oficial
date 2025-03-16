#include "TibiaImplantMatch.hpp"
#include <fstream>
#include "ImplantsException.hpp"
#include "ConvexHull.hpp"
#include "vtkCutter.h"
#include "vtkPlane.h"
#include "ImplantTools.hpp"
#include "vtkPlaneCollection.h"
#include "vtkClipClosedSurface.h"
#include "vtkExtractPolyDataGeometry.h"
#include "ImplantTools.hpp"
using namespace UKA::IMPLANTS;


inline double squareDistance2d(cv::Point2d P0, cv::Point2d P1)
{
	cv::Point2d diff = P1 - P0;
	return diff.dot(diff);
}

inline cv::Mat Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer transform)
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

inline cv::Mat Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer transform)
{
	itk::Vector< double, 3 > translate = transform->GetOffset();
	cv::Mat result(3, 1, CV_64FC1);

	result.at <double>(0, 0) = translate[0];
	result.at <double>(1, 0) = translate[1];
	result.at <double>(2, 0) = translate[2];
	return result;
}

inline cv::Mat GetRotateTransformAxisZ(Plane myPlane)
{
	myPlane.normalizeNormalVector();
	Point normalXY(0.0, 0.0, 1.0);
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

inline cv::Mat getTransformToRobot(Plane currentPlane, Point sagitalVector)
{
	cv::Mat rotateZ = GetRotateTransformAxisZ(currentPlane);
	cv::Mat rotateX = GetRotateTransformAxisX(sagitalVector, rotateZ);
	cv::Mat finalMatrix = rotateX * rotateZ;
	return finalMatrix;
}

TibiaImplantMatch::TibiaImplantMatch()
{
	isInit = false;
}

void TibiaImplantMatch::init(const TibiaImplant& implant, const Knee& knee)
{
	if (isInit == true)
	{
		throw ImplantExceptionCode::ALREADY_INITIALIZED_TIBIA_IMPLANT_MATCH;
	}
	this->implant = implant;
	this->knee = knee;
	//this->knee.setTibiaSlope(implant.getImplantInfo().slope);
	makeRotationMatrix();
	makeTranslationMatrix();
	isInit = true;
}

void TibiaImplantMatch::makeRotationMatrix()
{
	std::vector<cv::Point3d> implantVectors;
	std::vector<cv::Point3d> kneeVectors;

	Point implantNormal = implant.getTibiaNormalVector();
	Point implantAP = implant.getTibiaVectorAP();
	Point implantCross = implantNormal.cross(implantAP);
	implantVectors.push_back(implantNormal.ToCVPoint());
	implantVectors.push_back(implantAP.ToCVPoint());
	implantVectors.push_back(implantCross.ToCVPoint());

	Point kneeTibiaAxis = knee.getNormalVectorTibiaPlane();
	Point kneeTibiaAP = knee.getTibiaDirectVectorAP();
	Point kneeCross = kneeTibiaAxis.cross(kneeTibiaAP);

	kneeVectors.push_back(kneeTibiaAxis.ToCVPoint());
	kneeVectors.push_back(kneeTibiaAP.ToCVPoint());
	kneeVectors.push_back(kneeCross.ToCVPoint());

	cv::Mat implantMatrix = cv::Mat(implantVectors.size(), 3, CV_64F, implantVectors.data());
	cv::Mat kneeMatrix = cv::Mat(kneeVectors.size(), 3, CV_64F, kneeVectors.data());
	cv::Mat inverse = (implantMatrix.t()).inv();
	rotationMatrix = (kneeMatrix.t()) * inverse;
}

void TibiaImplantMatch::makeTranslationMatrix()
{
	double refDistance = 1000;
	Plane equisPlane = knee.getEquisPlaneTibia();
	Point vectorToSurgicalTEA = knee.getTibiaVectorToSurgicalSideTEA();
	
	Point vectorToSurgicalTEAProj = equisPlane.getProjectionVector(vectorToSurgicalTEA);
	vectorToSurgicalTEAProj.normalice();

	Point tibiaPlateau = knee.getMovePlateau(implant.getImplantInfo());
	//Point tibiaPlateau = knee.getTibiaCenterPointApAutomatic();
	translationMatrix = tibiaPlateau.ToMatPoint() - (rotationMatrix * implant.getPlateauRefPointDown().ToMatPoint());

	cv::Mat fullExtreme = rotationMatrix * implant.getExtremeSidePoint().ToMatPoint() + translationMatrix;
	Point currentPoint = equisPlane.getProjectionPoint(Point(fullExtreme));
	Point farPoint = currentPoint + refDistance * vectorToSurgicalTEAProj;

	vtkNew<vtkImplicitPolyDataDistance> polyDistance;
	polyDistance->SetInput(knee.GetTibiaPoly());

	Point borderPoint;

	double error = ImplantTools::GetInterceptionWithLine(polyDistance, currentPoint, farPoint, borderPoint);

	if (error <= 0.1)
	{
		double distance = ImplantTools::getDistanceBetweenPoints(borderPoint, farPoint);
		double moveDistance = refDistance - distance;
		//std::cout << "Error: " << error << " Move distance: " << moveDistance << std::endl;
		translationMatrix = translationMatrix + moveDistance * vectorToSurgicalTEA.ToMatPoint();
	}
}

Plane TibiaImplantMatch::transformPlane(const Plane& plane) const
{
	cv::Mat transformNormalVector = rotationMatrix * plane.getNormalVectorMat();
	cv::Mat transformPoint = (rotationMatrix * plane.getPointMat()) + translationMatrix;
	Plane transformPlane;
	transformPlane.init(Point(transformNormalVector), Point(transformPoint));
	return transformPlane;
}

Plane TibiaImplantMatch::finalTransformPlane(const Plane& plane, const itk::Rigid3DTransform<>::Pointer pTransform) const
{
	cv::Mat rotation = Rigid3DTransformToCVRotation(pTransform);
	cv::Mat translation = Rigid3DTransformToCVTranslation(pTransform);

	cv::Mat transformNormalVector = rotation * plane.getNormalVectorMat();
	cv::Mat transformPoint = (rotation * plane.getPointMat()) + translation;
	Plane transformPlane;
	transformPlane.init(Point(transformNormalVector), Point(transformPoint));
	return transformPlane;
}

Point TibiaImplantMatch::finalTransformPoint(const Point& pPoint, const itk::Rigid3DTransform<>::Pointer pTransform) const
{
	cv::Mat rotation = Rigid3DTransformToCVRotation(pTransform);
	cv::Mat translation = Rigid3DTransformToCVTranslation(pTransform);
	cv::Mat transformPoint = (rotation * pPoint.ToMatPoint()) + translation;
	return Point(transformPoint);
}

Plane TibiaImplantMatch::getTibiaPlane() const
{
	return transformPlane(implant.getTibiaPlane());
}

Point TibiaImplantMatch::transformImplantPoint(const Point& pPoint) const
{
	cv::Mat mat(3, 1, CV_64FC1);
	mat.at <double>(0, 0) = pPoint.x;
	mat.at <double>(1, 0) = pPoint.y;
	mat.at <double>(2, 0) = pPoint.z;

	cv::Mat transformPoint = (rotationMatrix * mat) + translationMatrix;
	return Point(transformPoint);
}

Knee TibiaImplantMatch::getKnee() const
{
	return knee;
}

std::vector<Point> TibiaImplantMatch::removeOutLiers(const std::vector<Point>& points, const Point& initPoint, const Point& lastPoint) const
{
	std::vector<Point> topBorder;
	std::list<Point> topBorderTemp;

	topBorder.push_back(initPoint);
	topBorderTemp.push_back(lastPoint);

	auto it1 = points.begin();
	auto it2 = points.end();

	for (; it1 != it2; ++it1)
	{
		topBorderTemp.push_back(*it1);
	}

	Point nearPoint, changePoint;
	changePoint = initPoint;
	do
	{
		if (deletePointsInsideRadius(topBorderTemp, changePoint, lastPoint, nearPoint))
		{
			changePoint = nearPoint;
			if (changePoint != initPoint)
			{
				topBorder.push_back(changePoint);
			}
		}

	} while (lastPoint != changePoint && topBorderTemp.size() > 0);

	return topBorder;
}

bool TibiaImplantMatch::deletePointsInsideRadius(std::list<Point>& points, const Point& centerPoint, const Point& diffPoint, Point& nearPoint, double radius) const
{
	bool result = false;
	double distance;
	double closeDistance = 9999999.0;
	std::list<Point>::iterator nearPointIt;
	auto it = points.begin();
	while (it != points.end())
	{
		distance = Line::getDistanceBetweenPoints(centerPoint, *it, true);
		if (distance < radius)
		{
			if ((*it) == diffPoint)
			{
				if (distance < closeDistance)
				{
					closeDistance = distance;
					nearPointIt = it;
					nearPoint = *it;
					result = true;
				}
				++it;
			}
			else
			{
				it = points.erase(it);
			}
		}
		else
		{
			if (distance < closeDistance)
			{
				closeDistance = distance;
				nearPointIt = it;
				nearPoint = *it;
				result = true;
			}
			++it;
		}
	}
	if (result == true)
	{
		points.erase(nearPointIt);
	}
	return result;
}

itk::Matrix< double, 3, 3 > TibiaImplantMatch::GetRotationMatrix() const
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

itk::Vector< double, 3 > TibiaImplantMatch::GetTranslationMatrix() const
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

vtkSmartPointer<vtkPolyData> TibiaImplantMatch::GetCuttingTibia() const
{
	Plane PlaneA = transformPlane(implant.getTibiaPlane());
	PlaneA.reverseByPoint(knee.getAnkleCenter());

	cv::Point3d planeNormal, planePoint;

	vtkNew<vtkPlane> vtkPlane;

	planeNormal = PlaneA.getNormalVector();
	planePoint = PlaneA.getPoint();
	vtkPlane->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
	vtkPlane->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);

	vtkNew<vtkPlaneCollection> tibiaPlanes;
	tibiaPlanes->AddItem(vtkPlane);

	vtkNew<vtkClipClosedSurface> tibiaClipper;
	tibiaClipper->SetInputData(knee.GetTibiaPoly());
	tibiaClipper->SetClippingPlanes(tibiaPlanes);
	tibiaClipper->Update();

	return tibiaClipper->GetOutput();
}

/*
std::vector<PointTypeITK> TibiaImplantMatch::GetHullPoints(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut, double distance, double distancePcl, double closeCurveLateral, double closeCurveMedial, int amount) const
{
	std::vector<PointTypeITK> hull;
	Plane myPlane = finalTransformPlane(implant.getTibiaPlane(), pTransformIn);

	Point normalTemp = knee.getAnkleCenter();
	Point myNormalFinal = myPlane.getProjectionPoint(normalTemp) - normalTemp;
	myPlane.fixNormalVector(myNormalFinal);

	double increaseVector = 100000.0;
	Point tubercle = myPlane.getProjectionPoint(knee.getTibiaTubercle());
	Point pcl = myPlane.getProjectionPoint(knee.getPclCenterPoint());
	Point tibiaCenter = myPlane.getProjectionPoint(knee.getTibiaCenterPointOnImplantAP(implant.getImplantInfo()));
	Point latPlateau = myPlane.getProjectionPoint(knee.getLateralPlateau());
	Point medPlateau = myPlane.getProjectionPoint(knee.getMedialPlateau());

	Point directVector = myPlane.getNormalVector();
	Point vectorAP = tubercle - pcl;
	vectorAP.normalice();
	Point vectorRobotAP = (-1.0) * vectorAP;

	cv::Mat myRotation = getTransformToRobot(myPlane, vectorRobotAP);

	Point vectorTrans;

	if (knee.getIsRight() == true)
	{
		vectorTrans = vectorAP.cross(directVector);
	}
	else
	{
		vectorTrans = directVector.cross(vectorAP);
	}

	vectorTrans.normalice();

	Point newPcl = pcl - increaseVector * vectorAP;

	Plane sagitalPlane, obliqueLatPlaneUp, obliqueMedPlaneUp, obliqueLatPlaneDown, obliqueMedPlaneDown;
	sagitalPlane.init(vectorTrans, pcl);

	if (closeCurveLateral > 1)
	{
		closeCurveLateral = 1.0;
	}

	if (closeCurveLateral < 0.05)
	{
		closeCurveLateral = 0.05;
	}

	if (closeCurveMedial > 1)
	{
		closeCurveMedial = 1.0;
	}

	if (closeCurveMedial < 0.05)
	{
		closeCurveMedial = 0.05;
	}

	Point obliquePointLat = tibiaCenter + vectorAP + 2.0 * closeCurveLateral * vectorTrans;
	//Point obliquePointMed = tibiaCenter + vectorAP - 3.0 * vectorTrans;

	Point refvectorMed = (tibiaCenter - vectorTrans) - (tibiaCenter + vectorAP);
	Point obliquePointMed = tibiaCenter + vectorAP + (closeCurveMedial + 0.5) * refvectorMed;
	obliqueLatPlaneDown = myPlane.getPerpendicularPlane(tibiaCenter, obliquePointLat);
	obliqueMedPlaneDown = myPlane.getPerpendicularPlane(tibiaCenter, obliquePointMed);

	if (obliqueLatPlaneDown.eval(obliquePointMed) < 0)
	{
		obliqueLatPlaneDown.reverse();
	}

	if (obliqueMedPlaneDown.eval(obliquePointLat) < 0)
	{
		obliqueMedPlaneDown.reverse();
	}

	Line lineTopContour(vectorTrans, newPcl);

	///////////////////////////////////////////////////////////////////

	vtkSmartPointer<vtkPolyData> contourMax = ImplantTools::getContours(knee.GetTibiaPoly(), myPlane.getNormalVector(), myPlane.getPoint());
	vtkSmartPointer<vtkPoints> vtkMyPoints = contourMax->GetPoints();
	int tVtkPointSize = vtkMyPoints->GetNumberOfPoints();


	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//ImplantTools::show(contourMax, knee.GetTibiaPoly());
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<Point> contourPoints;

	for (int i = 0; i < tVtkPointSize; i++)
	{
		double pnt[3];
		vtkMyPoints->GetPoint(i, pnt);
		Point currentPoint(pnt[0], pnt[1], pnt[2]);
		contourPoints.push_back(currentPoint);
	}

	if (contourPoints.size() < 20)
	{
		throw ImplantExceptionCode::NOT_ENOUGH_POINTS_ON_CUT_PLANE;
	}

	ConvexHull hullObj(contourPoints, myRotation);
	std::vector<Point> hullPoints = hullObj.GetConvexHull();

	if (hullPoints.size() < 4)
	{
		throw ImplantExceptionCode::CAN_NOT_DETERMINE_CONVEX_HULL_ON_TIBIA_CUT_PLANE;
	}

	int posLat, posMed;
	int tHullSize = hullPoints.size();
	Point hullCenter;

	hullCenter = ImplantTools::getHighestPointsOnTibia(hullPoints, sagitalPlane, lineTopContour, posLat, posMed);

	if (posLat == -1 || posMed == -1)
	{
		throw ImplantExceptionCode::CAN_NOT_DETERMINE_HIGHEST_POINTS_ON_BOTH_SIDE_OF_PCL;
	}

	int minPos, maxPos;

	if (posLat < posMed)
	{
		minPos = posLat;
		maxPos = posMed;
	}
	else
	{
		minPos = posMed;
		maxPos = posLat;
	}

	obliqueLatPlaneUp = myPlane.getPerpendicularPlane(hullCenter, hullPoints[posLat]);
	obliqueMedPlaneUp = myPlane.getPerpendicularPlane(hullCenter, hullPoints[posMed]);

	if (obliqueLatPlaneUp.eval(hullPoints[posMed]) < 0)
	{
		obliqueLatPlaneUp.reverse();
	}

	if (obliqueMedPlaneUp.eval(hullPoints[posLat]) < 0)
	{
		obliqueMedPlaneUp.reverse();
	}

	Point beginTop, lastTop;

	if (maxPos - minPos == 1)
	{
		beginTop = hullPoints[maxPos];
		lastTop = hullPoints[minPos];
		std::rotate(hullPoints.begin(), hullPoints.begin() + maxPos, hullPoints.end());
	}
	else
	{
		if (obliqueLatPlaneUp.eval(hullPoints[minPos + 1]) > 0 && obliqueMedPlaneUp.eval(hullPoints[minPos + 1]) > 0)
		{
			beginTop = hullPoints[maxPos];
			lastTop = hullPoints[minPos];
			std::rotate(hullPoints.begin(), hullPoints.begin() + maxPos, hullPoints.end());
		}
		else
		{
			beginTop = hullPoints[minPos];
			lastTop = hullPoints[maxPos];
			std::rotate(hullPoints.begin(), hullPoints.begin() + minPos, hullPoints.end());
		}
	}

	std::vector<Point> pclPoints;

	for (int i = 0; i < tVtkPointSize; i++)
	{
		double pnt[3];
		vtkMyPoints->GetPoint(i, pnt);
		Point currentPoint(pnt[0], pnt[1], pnt[2]);

		if (obliqueLatPlaneUp.eval(currentPoint) >= 0 && obliqueMedPlaneUp.eval(currentPoint) >= 0)
		{
			pclPoints.push_back(currentPoint);
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////
	//ImplantTools::show(contourMax, pclPoints);
	//ImplantTools::show(contourMax, contourMax);

	if (pclPoints.size() < 10)
	{
		throw ImplantExceptionCode::NOT_ENOUGH_POINTS_ON_PCL_BORDER;
	}

	std::vector<Point> hullConcave, finalHull;
	bool pclArea = false;
	for (int i = 0; i < tHullSize && pclArea == false; i++)
	{
		if (i < 2)
		{
			hullConcave.push_back(hullPoints[i]);
		}
		else if (!(obliqueLatPlaneUp.eval(hullPoints[i]) >= 0 && obliqueMedPlaneUp.eval(hullPoints[i]) >= 0))
		{
			hullConcave.push_back(hullPoints[i]);
		}
		else
		{
			pclArea = true;
		}
	}

	Plane topPlane = myPlane.getPerpendicularPlane(lastTop, beginTop);
	topPlane.reverseByPoint(hullCenter, false);

	Point xVector = myPlane.getNormalVector().cross(vectorRobotAP);
	xVector.normalice();
	xVector = -xVector;
	cv::Mat rotationPoly = getTransformToRobot(myPlane, xVector);

	//std::cout << "*********************************************************************" << std::endl;
	//ImplantTools::Poly tPoly = ImplantTools::polyFit(pclPoints, rotationPoly, 6, constraintObj);
	ImplantTools::Poly tPoly = ImplantTools::parabolaFitPCL(pclPoints, rotationPoly, beginTop, lastTop);
	//std::cout << "*********************************************************************" << std::endl;
	double maxX = tPoly.maxX;
	double minX = tPoly.minX;

	if (tPoly.isFine == false || maxX == minX)
	{
		throw ImplantExceptionCode::CAN_NOT_FIT_POLY_ON_PCL_BORDER;
	}

	std::vector<Point> hullConcaveTemp;
	double step = (maxX - minX) / 100;
	//double average = 20.;
	for (double i = minX; i < maxX + (step / 2.); i += step)
	{
		double y = tPoly.eval(i);
		Point temp(i, y, tPoly.Z);
		cv::Mat tempMat = rotationPoly.inv() * temp.ToMatPoint();
		hullConcaveTemp.push_back(Point(tempMat));
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//std::cout << tPoly.coeff[0] << ";  " << tPoly.coeff[1] << ";  " << tPoly.coeff[2] << std::endl;
	
	//auto vectorTemp1 = hullConcaveTemp;
	//vectorTemp1.push_back(beginTop);
	//vectorTemp1.push_back(lastTop);
	//ImplantTools::show(contourMax, pclPoints);
	//ImplantTools::show(contourMax, vectorTemp1);
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double reverseDistance = ImplantTools::getDistanceBetweenPoints(hullConcaveTemp[0], beginTop, true);
	double normalDistance = ImplantTools::getDistanceBetweenPoints(hullConcaveTemp[hullConcaveTemp.size() - 1], beginTop, true);

	if (reverseDistance < normalDistance)
	{
		std::reverse(hullConcaveTemp.begin(), hullConcaveTemp.end());
	}

	//////////////////////////////////////////////////////
	double slopeAngleLat = 90;
	double slopeAngleMed = 90;
	//double thicknessPlaneUp = 1.;
	//topPlane.movePlaneOnNormal(thicknessPlaneUp);
	Point oneTop = (lastTop + beginTop) / 2;
	sagitalPlane.movePlane(oneTop);

	Point sideLat = hullCenter + increaseVector * vectorTrans;
	Point sideMed = hullCenter - increaseVector * vectorTrans;

	Point interceptLat = topPlane.getInterceptionLinePoint(Line(topPlane.getNormalVector(), sideLat));
	Point interceptMed = topPlane.getInterceptionLinePoint(Line(topPlane.getNormalVector(), sideMed));

	int pDataSidePos;
	Line myLineObliqueLat = ImplantTools::GetSquareCornerFeatures(oneTop, interceptLat, sideLat, hullConcave, pDataSidePos, slopeAngleLat);
	interceptLat = topPlane.getInterceptionLinePoint(myLineObliqueLat);

	if (sagitalPlane.eval(beginTop) > 0)
	{
		while (pDataSidePos > 0)
		{
			hullConcave.erase(hullConcave.begin());
			pDataSidePos--;
		}

		Point lineVectorBegin = hullConcave[0];

		for (double i = 0.2; i < 0.9; i += 0.2)
		{
			Point temp = lineVectorBegin + i * (interceptLat - lineVectorBegin);
			hullConcave.insert(hullConcave.begin(), temp);
		}
	}
	else
	{
		while (hullConcave.size() > pDataSidePos + 1)
		{
			hullConcave.pop_back();
		}

		Point lineVectorBegin = hullConcave[pDataSidePos];

		for (double i = 0.2; i < 0.9; i += 0.2)
		{
			Point temp = lineVectorBegin + i * (interceptLat - lineVectorBegin);
			hullConcave.push_back(temp);
		}
	}

	Line myLineObliqueMed = ImplantTools::GetSquareCornerFeatures(oneTop, interceptMed, sideMed, hullConcave, pDataSidePos, slopeAngleMed);
	interceptMed = topPlane.getInterceptionLinePoint(myLineObliqueMed);

	if (sagitalPlane.eval(beginTop) < 0)
	{
		while (pDataSidePos > 0)
		{
			hullConcave.erase(hullConcave.begin());
			pDataSidePos--;
		}

		Point lineVectorBegin = hullConcave[0];

		for (double i = 0.2; i < 0.9; i += 0.2)
		{
			Point temp = lineVectorBegin + i * (interceptMed - lineVectorBegin);
			hullConcave.insert(hullConcave.begin(), temp);
		}
	}
	else
	{
		while (hullConcave.size() > pDataSidePos + 1)
		{
			hullConcave.pop_back();
		}

		Point lineVectorBegin = hullConcave[pDataSidePos];

		for (double i = 0.2; i < 0.9; i += 0.2)
		{
			Point temp = lineVectorBegin + i * (interceptMed - lineVectorBegin);
			hullConcave.push_back(temp);
		}
	}

	////////////////////////////////////////////////////

	int beginPos = -1, endPos = -1;

	beginPos = 0;
	endPos = hullConcaveTemp.size() - 1;

	
	//for (int i = 1; i < hullConcaveTemp.size(); i++)
	//{
	//	if (topPlane.eval(hullConcaveTemp[i]) < 0)
	//	{
	//		Point temp = topPlane.getInterceptionLinePoint(Line::makeLineWithPoints(hullConcaveTemp[i], hullConcaveTemp[i - 1]));
	//		hullConcaveTemp.insert(hullConcaveTemp.begin() + i, temp);
	//		beginPos = i;
	//		break;
	//	}
	//}

	//for (int i = hullConcaveTemp.size() - 2; i >= 0; i--)
	//{
	//	if (topPlane.eval(hullConcaveTemp[i]) < 0)
	//	{
	//		Point temp = topPlane.getInterceptionLinePoint(Line::makeLineWithPoints(hullConcaveTemp[i], hullConcaveTemp[i + 1]));
	//		hullConcaveTemp.insert(hullConcaveTemp.begin() + i + 1, temp);
	//		endPos = i + 1;
	//		break;
	//	}
	//}
	

	Point lineVectorBegin = hullConcave[hullConcave.size() - 1];
	Point lineVectorEnd = hullConcaveTemp[beginPos];

	int posBeginTopArea = hullConcave.size() - 1;

	for (double i = 0.2; i < 0.9; i += 0.2)
	{
		Point temp = lineVectorBegin + i * (lineVectorEnd - lineVectorBegin);
		hullConcave.push_back(temp);
	}

	int posBeginPCL = hullConcave.size();

	for (int i = beginPos; i <= endPos; i++)
	{
		hullConcave.push_back(hullConcaveTemp[i]);
	}

	int posEndPCL = hullConcave.size() - 1;

	lineVectorBegin = hullConcave[hullConcave.size() - 1];
	lineVectorEnd = hullConcave[0];

	for (double i = 0.2; i < 0.9; i += 0.2)
	{
		Point temp = lineVectorBegin + i * (lineVectorEnd - lineVectorBegin);
		hullConcave.push_back(temp);
	}

	
		//With this the last point ends in an upper corner and it be able to better move the curve top area.
	
	std::rotate(hullConcave.begin(), hullConcave.begin() + 1, hullConcave.end());
	int posEndTopArea = hullConcave.size() - 1;
	posBeginTopArea = posBeginTopArea - 1;
	posBeginPCL = posBeginPCL - 1;
	posEndPCL = posEndPCL - 1;

	Point computePoint;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//auto poly = ImplantTools::getPolyLine(hullConcave);
	//std::vector<Point> tExample = { hullConcave[posBeginPCL], hullConcave[posEndPCL], hullConcave[posEndTopArea], hullConcave[posBeginTopArea] };
	//ImplantTools::show(poly, tExample);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (int i = 0; i < hullConcave.size(); i++)
	{
		if (i >= posBeginTopArea && i <= posEndTopArea)
		{
			computePoint = hullConcave[i] + distancePcl * topPlane.getNormalVector();
			finalHull.push_back(computePoint);
		}
		else if (distance != 0)
		{
			if (i + 1 != hullConcave.size())
			{
				computePoint = movePointAtNormal(hullConcave[i], hullConcave[i + 1], myRotation, distance);
				finalHull.push_back(computePoint);
			}
		}
		else
		{
			finalHull.push_back(hullConcave[i]);
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//auto poly2 = ImplantTools::getPolyLine(hullConcave);
	//ImplantTools::show(poly2, finalHull, true);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<Point> concaveSpline = ConvexHull::interpolateSpline(finalHull, amount);
	int initPos = -1;

	for (int i = 0; i < concaveSpline.size(); i++)
	{
		if (obliqueLatPlaneDown.eval(concaveSpline[i]) >= 0 && obliqueMedPlaneDown.eval(concaveSpline[i]) >= 0)
		{
			initPos = i;
			break;
		}
	}

	if (initPos == -1)
	{
		initPos = 0;
	}

	if (initPos > 0)
	{
		std::rotate(concaveSpline.begin(), concaveSpline.begin() + initPos, concaveSpline.end());
	}

	std::vector<Point> finalSplinePoints;

	for (int i = 0; i < concaveSpline.size(); i++)
	{
		if (!(obliqueLatPlaneDown.eval(concaveSpline[i]) >= 0 && obliqueMedPlaneDown.eval(concaveSpline[i]) >= 0))
		{
			finalSplinePoints.push_back(concaveSpline[i]);
		}
	}

	Point initExtreme, endExtreme;

	if (finalSplinePoints.size() > 1)
	{
		initExtreme = finalSplinePoints[0];
		endExtreme = finalSplinePoints[finalSplinePoints.size() - 1];

		cv::Mat initExtremeMat, endExtremeMat;

		initExtremeMat = myRotation * initExtreme.ToMatPoint();
		endExtremeMat = myRotation * endExtreme.ToMatPoint();

		initExtreme = Point(initExtremeMat);
		endExtreme = Point(endExtremeMat);

		if (endExtreme.y < initExtreme.y)
		{
			std::reverse(finalSplinePoints.begin(), finalSplinePoints.end());
		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//ImplantTools::show(contourMax, finalSplinePoints);

	hull = increaseVectorToAmount(finalSplinePoints, amount);

	itk::Vector< double, 3 > translate;

	if (tHullSize > 0)
	{
		hullCenter = sagitalPlane.getProjectionPoint(hullCenter);
		hullCenter = myPlane.getProjectionPoint(hullCenter);
		Point tTemp = myRotation * (hullCenter.ToMatPoint());
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
*/

/*
std::vector<PointTypeITK> TibiaImplantMatch::GetHullPointsOld(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut, double distance, double distancePcl, int amount) const
{
	std::vector<PointTypeITK> hull;
	Plane myPlane = finalTransformPlane(implant.getTibiaPlane(), pTransformIn);

	Point normalTemp = knee.getAnkleCenter();
	Point myNormalFinal = myPlane.getProjectionPoint(normalTemp) - normalTemp;
	myPlane.fixNormalVector(myNormalFinal);

	std::vector<Point> tibiaProjected;
	double increaseVector = 100000.0;
	Point tubercle = myPlane.getProjectionPoint(knee.getTibiaTubercle());
	Point pcl = myPlane.getProjectionPoint(knee.getPclCenterPoint());
	Point tibiaCenter = myPlane.getProjectionPoint(knee.getTibiaCenterPointOnImplantAP(implant.getImplantInfo()));
	Point latPlateau = myPlane.getProjectionPoint(knee.getLateralPlateau());
	Point medPlateau = myPlane.getProjectionPoint(knee.getMedialPlateau());
	Point plateau;

	Point directVector = myPlane.getNormalVector();
	Point vectorAP = tubercle - pcl;
	vectorAP.normalice();
	Point vectorRobotAP = (-1.0) * vectorAP;

	cv::Mat myRotation = getTransformToRobot(myPlane, vectorRobotAP);

	Point vectorTrans = directVector.cross(vectorAP);

	Point newTubercle = tubercle + increaseVector * vectorAP;
	Point newPcl = pcl - increaseVector * vectorAP;

	Plane sagitalPlane, obliqueLatPlaneUp, obliqueMedPlaneUp, obliqueLatPlaneDown, obliqueMedPlaneDown;
	sagitalPlane.init(vectorTrans, pcl);
	//sagitalPlane.normalizeNormalVector();

	if (sagitalPlane.isPointBelongToPlane(latPlateau) == true)
	{
		if (sagitalPlane.eval(medPlateau) > 0)
		{
			sagitalPlane.reverse();
		}
		vectorTrans = sagitalPlane.getProjectionPoint(medPlateau) - medPlateau;
	}
	else
	{
		if (sagitalPlane.eval(latPlateau) < 0)
		{
			sagitalPlane.reverse();
		}
		vectorTrans = latPlateau - sagitalPlane.getProjectionPoint(latPlateau);
	}

	vectorTrans.normalice();

	Point obliquePointLat = tibiaCenter + vectorAP + 1.2 * vectorTrans;
	Point obliquePointMed = tibiaCenter + vectorAP - 3.0 * vectorTrans;
	obliqueLatPlaneDown = myPlane.getPerpendicularPlane(tibiaCenter, obliquePointLat);
	obliqueMedPlaneDown = myPlane.getPerpendicularPlane(tibiaCenter, obliquePointMed);

	if (obliqueLatPlaneDown.eval(obliquePointMed) < 0)
	{
		obliqueLatPlaneDown.reverse();
	}

	if (obliqueMedPlaneDown.eval(obliquePointLat) < 0)
	{
		obliqueMedPlaneDown.reverse();
	}

	Line lineTopLat(vectorTrans, newPcl);
	Line lineTopMed(vectorTrans, newPcl);

	double topLateralDistance = -1.0;
	double topMedialDistance = -1.0;
	double topLateral, topMedial;
	Point topLateralPoint, topMedialPoint;

	vtkSmartPointer<vtkPolyData> contour = getContour(knee.GetTibiaPoly(), myPlane.getNormalVector(), myPlane.getPoint());
	vtkSmartPointer<vtkPoints> vtkMyPoints = contour->GetPoints();
	int tVtkPointSize = vtkMyPoints->GetNumberOfPoints();

	Point tempPoint;
	double myEpsilon = 0.01;
	Point midPlanePoint;
	int cont = 0;

	for (int i = 0; i < tVtkPointSize; i++)
	{
		double pnt[3];
		vtkMyPoints->GetPoint(i, pnt);
		tempPoint = Point(pnt[0], pnt[1], pnt[2]);

		midPlanePoint = midPlanePoint + tempPoint;
		cont++;
		if (obliqueLatPlaneDown.eval(tempPoint) > 0 && obliqueMedPlaneDown.eval(tempPoint) > 0)
		{
			continue;
		}

		tibiaProjected.push_back(tempPoint);

		if (sagitalPlane.eval(tempPoint) > 0)
		{
			topLateral = 1.0 / (lineTopLat.getDistanceFromPoint(tempPoint) + myEpsilon);
			if (topLateral > topLateralDistance)
			{
				topLateralDistance = topLateral;
				topLateralPoint = tempPoint;
			}
		}
		else
		{
			topMedial = 1.0 / (lineTopMed.getDistanceFromPoint(tempPoint) + myEpsilon);
			if (topMedial > topMedialDistance)
			{
				topMedialDistance = topMedial;
				topMedialPoint = tempPoint;
			}
		}
	}
	if (tibiaProjected.size() < 15)
	{
		throw ImplantsException("There are not enough points on the cutting plane.");
	}

	midPlanePoint = midPlanePoint / double(cont);

	std::vector<Point>::iterator it1 = tibiaProjected.begin();
	std::vector<Point>::iterator it2 = tibiaProjected.end();

	ConvexHull hullObj(tibiaProjected, myRotation);
	std::vector<Point> hullPoints = hullObj.GetConvexHull();

	if (hullPoints.size() == 0)
	{
		throw ImplantsException("The convex hull on tibia cut plane could not be determined.");
	}

	Point hullLateral, hullMedial;
	Point hullCenter;
	double distanceHullLat, distanceHullMed;
	distanceHullLat = 1000000.0;
	distanceHullMed = 1000000.0;
	for (int i = 0; i < hullPoints.size(); i++)
	{
		hullCenter = hullCenter + hullPoints[i];
		double distLat = Line::getDistanceBetweenPoints(topLateralPoint, hullPoints[i], true);
		double distMed = Line::getDistanceBetweenPoints(topMedialPoint, hullPoints[i], true);

		if (distLat < distanceHullLat)
		{
			distanceHullLat = distLat;
			hullLateral = hullPoints[i];
		}

		if (distMed < distanceHullMed)
		{
			distanceHullMed = distMed;
			hullMedial = hullPoints[i];
		}
	}

	hullCenter = hullCenter / double(hullPoints.size());

	std::vector<Point> listMainPoints = { hullLateral, hullMedial, hullCenter };

	obliqueLatPlaneUp = myPlane.getPerpendicularPlane(listMainPoints[2], listMainPoints[0]);
	obliqueMedPlaneUp = myPlane.getPerpendicularPlane(listMainPoints[2], listMainPoints[1]);

	if (obliqueLatPlaneUp.eval(listMainPoints[1]) < 0)
	{
		obliqueLatPlaneUp.reverse();
	}

	if (obliqueMedPlaneUp.eval(listMainPoints[0]) < 0)
	{
		obliqueMedPlaneUp.reverse();
	}

	int initPos = -1;

	for (int i = 0; i < hullPoints.size(); i++)
	{
		if (obliqueLatPlaneUp.eval(hullPoints[i]) >= 0 && obliqueMedPlaneUp.eval(hullPoints[i]) >= 0)
		{
			initPos = i;
			break;
		}
	}
	if (initPos == -1)
	{
		initPos = 0;
	}

	if (initPos > 0)
	{
		std::rotate(hullPoints.begin(), hullPoints.begin() + initPos, hullPoints.end());
	}

	std::vector<Point> concavePoints;
	it1 = tibiaProjected.begin();
	it2 = tibiaProjected.end();

	for (; it1 != it2; ++it1)
	{
		if (obliqueLatPlaneUp.eval(*it1) >= 0 && obliqueMedPlaneUp.eval(*it1) >= 0)
		{
			concavePoints.push_back((*it1));
		}
	}

	std::vector<Point> concaveHull, concaveHullTemp, finalConcaveTop, ConcaveTop;

	concaveHullTemp = removeOutLiers(concavePoints, listMainPoints[0], listMainPoints[1]);

	if (concaveHullTemp.size() == 0)
	{
		throw ImplantsException("Outliers points could not be removed in top border.");
	}

	Point computePoint;

	for (int i = 0; i < hullPoints.size() - 1; i++)
	{
		if (!(obliqueLatPlaneUp.eval(hullPoints[i]) >= 0 && obliqueMedPlaneUp.eval(hullPoints[i]) >= 0))
		{
			computePoint = movePointAtNormal(hullPoints[i], hullPoints[i + 1], myRotation, distance);
			concaveHull.push_back(computePoint);
		}

		if (i == hullPoints.size() - 2)
		{
			if (!(obliqueLatPlaneUp.eval(hullPoints[i + 1]) >= 0 && obliqueMedPlaneUp.eval(hullPoints[i + 1]) >= 0))
			{
				computePoint = movePointAtNormal(hullPoints[i], hullPoints[i + 1], myRotation, distance, true);
				concaveHull.push_back(computePoint);
			}
		}
	}

	for (int i = 0; i < concaveHullTemp.size(); i += 5)
	{
		ConcaveTop.push_back(concaveHullTemp[i]);
	}

	if (concaveHullTemp.size() % 5 == 0 || concaveHullTemp.size() % 5 > 1)
	{
		ConcaveTop.push_back(concaveHullTemp[concaveHullTemp.size() - 1]);
	}

	finalConcaveTop = ConvexHull::ReduceConvexHull(ConcaveTop, 100);

	double initDist = Line::getDistanceBetweenPoints(finalConcaveTop[0], concaveHull[concaveHull.size() - 1], true);
	double lastDist = Line::getDistanceBetweenPoints(finalConcaveTop[finalConcaveTop.size() - 1], concaveHull[concaveHull.size() - 1], true);

	if (initDist < lastDist)
	{
		for (int i = 0; i < finalConcaveTop.size() - 1; i++)
		{

			computePoint = movePointAtNormal(finalConcaveTop[i], finalConcaveTop[i + 1], myRotation, distancePcl);

			concaveHull.push_back(computePoint);

			if (i == finalConcaveTop.size() - 2)
			{
				computePoint = movePointAtNormal(finalConcaveTop[i], finalConcaveTop[i + 1], myRotation, distancePcl, true);
				concaveHull.push_back(computePoint);
			}
		}
	}
	else
	{
		for (int i = finalConcaveTop.size() - 1; i > 0; i--)
		{
			computePoint = movePointAtNormal(finalConcaveTop[i], finalConcaveTop[i - 1], myRotation, distancePcl);

			concaveHull.push_back(computePoint);

			if (i == 1)
			{
				computePoint = movePointAtNormal(finalConcaveTop[i], finalConcaveTop[i - 1], myRotation, distancePcl, true);
				concaveHull.push_back(computePoint);
			}
		}
	}

	std::vector<Point> concaveSpline = ConvexHull::interpolateSpline(concaveHull, amount);
	initPos = -1;

	for (int i = 0; i < concaveSpline.size(); i++)
	{
		if (obliqueLatPlaneDown.eval(concaveSpline[i]) >= 0 && obliqueMedPlaneDown.eval(concaveSpline[i]) >= 0)
		{
			initPos = i;
			break;
		}
	}

	if (initPos == -1)
	{
		initPos = 0;
	}

	if (initPos > 0)
	{
		std::rotate(concaveSpline.begin(), concaveSpline.begin() + initPos, concaveSpline.end());
	}

	std::vector<Point> finalSplinePoints;

	for (int i = 0; i < concaveSpline.size(); i++)
	{
		if (!(obliqueLatPlaneDown.eval(concaveSpline[i]) >= 0 && obliqueMedPlaneDown.eval(concaveSpline[i]) >= 0))
		{
			finalSplinePoints.push_back(concaveSpline[i]);
		}
	}

	Point initExtreme, endExtreme;

	if (finalSplinePoints.size() > 1)
	{
		initExtreme = finalSplinePoints[0];
		endExtreme = finalSplinePoints[finalSplinePoints.size() - 1];

		cv::Mat initExtremeMat, endExtremeMat;

		initExtremeMat = myRotation * initExtreme.ToMatPoint();
		endExtremeMat = myRotation * endExtreme.ToMatPoint();

		initExtreme = Point(initExtremeMat);
		endExtreme = Point(endExtremeMat);

		if (endExtreme.y < initExtreme.y)
		{
			std::reverse(finalSplinePoints.begin(), finalSplinePoints.end());
		}
	}

	hull = increaseVectorToAmount(finalSplinePoints, amount);

	itk::Vector< double, 3 > translate;
	if (tibiaProjected.size() > 0)
	{
		midPlanePoint = sagitalPlane.getProjectionPoint(midPlanePoint);
		midPlanePoint = myPlane.getProjectionPoint(midPlanePoint);
		Point tTemp = myRotation * (midPlanePoint.ToMatPoint());
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

	//std::ofstream outfile("tibia_plane.txt");
	//double a, b, c;
	//Point example1, P1, P2;

	//for (int i = 0; i < hull.size(); i++)
	//{
	//    if (i == hull.size() - 1)
	//    {
	//        P1 = printVector[i];
	//        P2 = printVector[0];
	//    }
	//    else
	//    {
	//        PointTypeITK aa, bb;
	//        aa = hull[i];
	//        bb = hull[i + 1];
	//        P1 = Point(aa[0], aa[1], aa[2]);
	//        P2 = Point(bb[0], bb[1], bb[2]);
	//    }

	//    for (double j = 0; j < 1001.0; j++)
	//    {
	//        example1 = (P1 + (j / 1000)*(P2 - P1));
	//        outfile << example1.x << " " << example1.y << " " << " " << example1.z << "\n";
	//    }

	//    PointTypeITK aa;
	//    aa = reducePointsSpline[i];
	//    P1 = Point(aa[0], aa[1], aa[2]);
	//    cv::Mat pp = myRotation * P1.ToMatPoint();
	//    example1 = Point(pp);
	//    outfile << example1.x << " " << example1.y << " " << " " << example1.z << "\n";
	//}
	//outfile.close();

	return hull;
}

*/

Point TibiaImplantMatch::movePointAtNormal(const Point& movePoint, const Point& nextPoint, const cv::Mat& rotationZ, double distance, bool changeMove, bool clockWise) const
{
	cv::Mat moveMat = rotationZ * movePoint.ToMatPoint();
	cv::Mat nextMat = rotationZ * nextPoint.ToMatPoint();

	Point moveP = Point(moveMat);
	Point nextP = Point(nextMat);
	Point vector;

	if (clockWise == false)
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

std::vector<PointTypeITK> TibiaImplantMatch::increaseVectorToAmount(const std::vector<Point>& points, int amount) const
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

	/*int lastPos = 0;
	int stillPoints = amount - points.size();
	int interAmount = stillPoints / interPoints;
	int rest = stillPoints % interPoints;
	Point a, b, c;
	double coef = 1.0 / double(interPoints + 1);
	for (int i = 0; i < points.size() - 1; i++)
	{
		lastPos = i;
		a = points[i];
		b = points[i + 1];
		result.push_back(a.ToITKPoint());
		for (int j = 1; j <= interPoints; j++)
		{
			c = a + double(j) * coef * (b - a);
			result.push_back(c.ToITKPoint());
		}
		interAmount--;
		if (interAmount == 0)
		{
			break;
		}
	}*/

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

vtkSmartPointer<vtkPolyData> TibiaImplantMatch::getContour(const vtkSmartPointer<vtkPolyData> poly, const Point& pNormal, const Point& pPoint) const
{
	vtkNew<vtkPlane> cutPlane;
	cutPlane->SetOrigin(pPoint.x, pPoint.y, pPoint.z);
	cutPlane->SetNormal(pNormal.x, pNormal.y, pNormal.z);

	vtkNew<vtkCutter> cutter;
	cutter->SetInputData(poly);
	cutter->SetCutFunction(cutPlane);
	cutter->Update();

	auto contour = cutter->GetOutput();

	return contour;
}

std::vector<PointTypeITK> TibiaImplantMatch::GetHullPoints(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut, double distance, double distancePcl, double distanceSide, double closeCurveLateral, double closeCurveMedial, int amount) const
{
	std::vector<PointTypeITK> hull;
	Plane myPlane = finalTransformPlane(implant.getTibiaPlane(), pTransformIn);

	Point normalTemp = knee.getAnkleCenter();
	Point myNormalFinal = myPlane.getProjectionPoint(normalTemp) - normalTemp;
	myPlane.fixNormalVector(myNormalFinal);

	double increaseVector = 100000.0;
	Point tubercle = myPlane.getProjectionPoint(knee.getTibiaTubercle());
	Point pcl = myPlane.getProjectionPoint(knee.getPclCenterPoint());
	Point tibiaCenter = myPlane.getProjectionPoint(knee.getTibiaCenterPointOnImplantAP(implant.getImplantInfo()));
	Point latPlateau = myPlane.getProjectionPoint(knee.getLateralPlateau());
	Point medPlateau = myPlane.getProjectionPoint(knee.getMedialPlateau());

	Point directVector = myPlane.getNormalVector();
	Point vectorAP = tubercle - pcl;
	vectorAP.normalice();
	Point vectorRobotAP = (-1.0) * vectorAP;

	cv::Mat myRotation = getTransformToRobot(myPlane, vectorRobotAP);

	Point vectorTrans;

	if (knee.getIsRight() == true)
	{
		vectorTrans = vectorAP.cross(directVector);
	}
	else
	{
		vectorTrans = directVector.cross(vectorAP);
	}

	vectorTrans.normalice();

	Point newPcl = pcl - increaseVector * vectorAP;

	Plane sagitalPlane, obliqueLatPlaneUp, obliqueMedPlaneUp, obliqueLatPlaneDown, obliqueMedPlaneDown;
	sagitalPlane.init(vectorTrans, pcl);

	if (closeCurveLateral > 1)
	{
		closeCurveLateral = 1.0;
	}

	if (closeCurveLateral < 0.05)
	{
		closeCurveLateral = 0.05;
	}

	if (closeCurveMedial > 1)
	{
		closeCurveMedial = 1.0;
	}

	if (closeCurveMedial < 0.05)
	{
		closeCurveMedial = 0.05;
	}

	Point obliquePointLat = tibiaCenter + vectorAP + 2.0 * closeCurveLateral * vectorTrans;
	//Point obliquePointMed = tibiaCenter + vectorAP - 3.0 * vectorTrans;

	Point refvectorMed = (tibiaCenter - vectorTrans) - (tibiaCenter + vectorAP);
	Point obliquePointMed = tibiaCenter + vectorAP + (closeCurveMedial + 0.5) * refvectorMed;
	obliqueLatPlaneDown = myPlane.getPerpendicularPlane(tibiaCenter, obliquePointLat);
	obliqueMedPlaneDown = myPlane.getPerpendicularPlane(tibiaCenter, obliquePointMed);

	if (obliqueLatPlaneDown.eval(obliquePointMed) < 0)
	{
		obliqueLatPlaneDown.reverse();
	}

	if (obliqueMedPlaneDown.eval(obliquePointLat) < 0)
	{
		obliqueMedPlaneDown.reverse();
	}

	Line lineTopContour(vectorTrans, newPcl);

	///////////////////////////////////////////////////////////////////

	vtkSmartPointer<vtkPolyData> contourMax = ImplantTools::getContours(knee.GetTibiaPoly(), myPlane.getNormalVector(), myPlane.getPoint());
	vtkSmartPointer<vtkPoints> vtkMyPoints = contourMax->GetPoints();
	int tVtkPointSize = vtkMyPoints->GetNumberOfPoints();


	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//ImplantTools::show(contourMax, knee.GetTibiaPoly());
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<Point> contourPoints;

	for (int i = 0; i < tVtkPointSize; i++)
	{
		double pnt[3];
		vtkMyPoints->GetPoint(i, pnt);
		Point currentPoint(pnt[0], pnt[1], pnt[2]);
		contourPoints.push_back(currentPoint);
	}

	if (contourPoints.size() < 20)
	{
		throw ImplantExceptionCode::NOT_ENOUGH_POINTS_ON_CUT_PLANE;
	}

	ConvexHull hullObj(contourPoints, myRotation);
	std::vector<Point> hullPoints = hullObj.GetConvexHull();
	if (hullPoints.size() < 4)
	{
		throw ImplantExceptionCode::CAN_NOT_DETERMINE_CONVEX_HULL_ON_TIBIA_CUT_PLANE;
	}

	int posLat, posMed;
	int tHullSize = hullPoints.size();
	Point hullCenter;

	hullCenter = ImplantTools::getHighestPointsOnTibia(hullPoints, sagitalPlane, lineTopContour, posLat, posMed);

	if (posLat == -1 || posMed == -1)
	{
		throw ImplantExceptionCode::CAN_NOT_DETERMINE_HIGHEST_POINTS_ON_BOTH_SIDE_OF_PCL;
	}

	int minPos, maxPos;

	if (posLat < posMed)
	{
		minPos = posLat;
		maxPos = posMed;
	}
	else
	{
		minPos = posMed;
		maxPos = posLat;
	}

	obliqueLatPlaneUp = myPlane.getPerpendicularPlane(hullCenter, hullPoints[posLat]);
	obliqueMedPlaneUp = myPlane.getPerpendicularPlane(hullCenter, hullPoints[posMed]);

	if (obliqueLatPlaneUp.eval(hullPoints[posMed]) < 0)
	{
		obliqueLatPlaneUp.reverse();
	}

	if (obliqueMedPlaneUp.eval(hullPoints[posLat]) < 0)
	{
		obliqueMedPlaneUp.reverse();
	}

	Point beginTop, lastTop;

	if (maxPos - minPos == 1)
	{
		beginTop = hullPoints[maxPos];
		lastTop = hullPoints[minPos];
		std::rotate(hullPoints.begin(), hullPoints.begin() + maxPos, hullPoints.end());
	}
	else
	{
		if (obliqueLatPlaneUp.eval(hullPoints[minPos + 1]) > 0 && obliqueMedPlaneUp.eval(hullPoints[minPos + 1]) > 0)
		{
			beginTop = hullPoints[maxPos];
			lastTop = hullPoints[minPos];
			std::rotate(hullPoints.begin(), hullPoints.begin() + maxPos, hullPoints.end());
		}
		else
		{
			beginTop = hullPoints[minPos];
			lastTop = hullPoints[maxPos];
			std::rotate(hullPoints.begin(), hullPoints.begin() + minPos, hullPoints.end());
		}
	}

	std::vector<Point> pclPoints;

	for (int i = 0; i < tVtkPointSize; i++)
	{
		double pnt[3];
		vtkMyPoints->GetPoint(i, pnt);
		Point currentPoint(pnt[0], pnt[1], pnt[2]);

		if (obliqueLatPlaneUp.eval(currentPoint) >= 0 && obliqueMedPlaneUp.eval(currentPoint) >= 0)
		{
			pclPoints.push_back(currentPoint);
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////
	//ImplantTools::show(contourMax, pclPoints);
	//ImplantTools::show(contourMax, contourMax);
	//////////////////////////////////////////////////////////////////////////////////

	if (pclPoints.size() < 10)
	{
		throw ImplantExceptionCode::NOT_ENOUGH_POINTS_ON_PCL_BORDER;
	}

	std::vector<Point> hullConcave, finalHull;
	bool pclArea = false;
	for (int i = 0; i < tHullSize && pclArea == false; i++)
	{
		if (i < 2)
		{
			hullConcave.push_back(hullPoints[i]);
		}
		else if (!(obliqueLatPlaneUp.eval(hullPoints[i]) >= 0 && obliqueMedPlaneUp.eval(hullPoints[i]) >= 0))
		{
			hullConcave.push_back(hullPoints[i]);
		}
		else
		{
			pclArea = true;
		}
	}

	Plane topPlane = myPlane.getPerpendicularPlane(lastTop, beginTop);
	topPlane.reverseByPoint(hullCenter, false);

	Point xVector = myPlane.getNormalVector().cross(vectorRobotAP);
	xVector.normalice();
	xVector = -xVector;
	cv::Mat rotationPoly = getTransformToRobot(myPlane, xVector);

	//std::cout << "*********************************************************************" << std::endl;
	//ImplantTools::Poly tPoly = ImplantTools::polyFit(pclPoints, rotationPoly, 6, constraintObj);
	ImplantTools::Poly tPoly = ImplantTools::parabolaFitPCL(pclPoints, rotationPoly, beginTop, lastTop);
	//std::cout << "*********************************************************************" << std::endl;
	double maxX = tPoly.maxX;
	double minX = tPoly.minX;

	if (tPoly.isFine == false || maxX == minX)
	{
		throw ImplantExceptionCode::CAN_NOT_FIT_POLY_ON_PCL_BORDER;
	}

	std::vector<Point> hullConcaveTemp;
	double step = (maxX - minX) / 100;
	//double average = 20.;
	for (double i = minX; i < maxX + (step / 2.); i += step)
	{
		double y = tPoly.eval(i);
		Point temp(i, y, tPoly.Z);
		cv::Mat tempMat = rotationPoly.inv() * temp.ToMatPoint();
		hullConcaveTemp.push_back(Point(tempMat));
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*std::cout << tPoly.coeff[0] << ";  " << tPoly.coeff[1] << ";  " << tPoly.coeff[2] << std::endl;*/
	/*
	auto vectorTemp1 = hullConcaveTemp;
	vectorTemp1.push_back(beginTop);
	vectorTemp1.push_back(lastTop);
	ImplantTools::show(contourMax, pclPoints);
	ImplantTools::show(contourMax, vectorTemp1);
	*/
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double reverseDistance = ImplantTools::getDistanceBetweenPoints(hullConcaveTemp[0], beginTop, true);
	double normalDistance = ImplantTools::getDistanceBetweenPoints(hullConcaveTemp[hullConcaveTemp.size() - 1], beginTop, true);

	if (reverseDistance < normalDistance)
	{
		std::reverse(hullConcaveTemp.begin(), hullConcaveTemp.end());
	}

	//////////////////////////////////////////////////////
	double slopeAngleLat = 90;
	double slopeAngleMed = 90;
	//double thicknessPlaneUp = 1.;
	//topPlane.movePlaneOnNormal(thicknessPlaneUp);
	Point oneTop = (lastTop + beginTop) / 2;
	sagitalPlane.movePlane(oneTop);

	Point sideLat = hullCenter + increaseVector * vectorTrans;
	Point sideMed = hullCenter - increaseVector * vectorTrans;

	Point interceptLat = topPlane.getInterceptionLinePoint(Line(topPlane.getNormalVector(), sideLat));
	Point interceptMed = topPlane.getInterceptionLinePoint(Line(topPlane.getNormalVector(), sideMed));

	int pDataSidePos;
	Line myLineObliqueLat = ImplantTools::GetSquareCornerFeatures(oneTop, interceptLat, sideLat, hullConcave, pDataSidePos, slopeAngleLat);
	interceptLat = topPlane.getInterceptionLinePoint(myLineObliqueLat);

	if (sagitalPlane.eval(beginTop) > 0)
	{
		while (pDataSidePos > 0)
		{
			hullConcave.erase(hullConcave.begin());
			pDataSidePos--;
		}

		Point lineVectorBegin = hullConcave[0];

		for (double i = 0.2; i < 0.9; i += 0.2)
		{
			Point temp = lineVectorBegin + i * (interceptLat - lineVectorBegin);
			hullConcave.insert(hullConcave.begin(), temp);
		}
	}
	else
	{
		while (hullConcave.size() > pDataSidePos + 1)
		{
			hullConcave.pop_back();
		}

		Point lineVectorBegin = hullConcave[pDataSidePos];

		for (double i = 0.2; i < 0.9; i += 0.2)
		{
			Point temp = lineVectorBegin + i * (interceptLat - lineVectorBegin);
			hullConcave.push_back(temp);
		}
	}

	Line myLineObliqueMed = ImplantTools::GetSquareCornerFeatures(oneTop, interceptMed, sideMed, hullConcave, pDataSidePos, slopeAngleMed);
	interceptMed = topPlane.getInterceptionLinePoint(myLineObliqueMed);

	if (sagitalPlane.eval(beginTop) < 0)
	{
		while (pDataSidePos > 0)
		{
			hullConcave.erase(hullConcave.begin());
			pDataSidePos--;
		}

		Point lineVectorBegin = hullConcave[0];

		for (double i = 0.2; i < 0.9; i += 0.2)
		{
			Point temp = lineVectorBegin + i * (interceptMed - lineVectorBegin);
			hullConcave.insert(hullConcave.begin(), temp);
		}
	}
	else
	{
		while (hullConcave.size() > pDataSidePos + 1)
		{
			hullConcave.pop_back();
		}

		Point lineVectorBegin = hullConcave[pDataSidePos];

		for (double i = 0.2; i < 0.9; i += 0.2)
		{
			Point temp = lineVectorBegin + i * (interceptMed - lineVectorBegin);
			hullConcave.push_back(temp);
		}
	}

	////////////////////////////////////////////////////

	int beginPos = -1, endPos = -1;

	beginPos = 0;
	endPos = hullConcaveTemp.size() - 1;

	/*
	for (int i = 1; i < hullConcaveTemp.size(); i++)
	{
		if (topPlane.eval(hullConcaveTemp[i]) < 0)
		{
			Point temp = topPlane.getInterceptionLinePoint(Line::makeLineWithPoints(hullConcaveTemp[i], hullConcaveTemp[i - 1]));
			hullConcaveTemp.insert(hullConcaveTemp.begin() + i, temp);
			beginPos = i;
			break;
		}
	}

	for (int i = hullConcaveTemp.size() - 2; i >= 0; i--)
	{
		if (topPlane.eval(hullConcaveTemp[i]) < 0)
		{
			Point temp = topPlane.getInterceptionLinePoint(Line::makeLineWithPoints(hullConcaveTemp[i], hullConcaveTemp[i + 1]));
			hullConcaveTemp.insert(hullConcaveTemp.begin() + i + 1, temp);
			endPos = i + 1;
			break;
		}
	}
	*/

	Point lineVectorBegin = hullConcave[hullConcave.size() - 1];
	Point lineVectorEnd = hullConcaveTemp[beginPos];

	int posBeginTopArea = hullConcave.size() - 1;

	for (double i = 0.2; i < 0.9; i += 0.2)
	{
		Point temp = lineVectorBegin + i * (lineVectorEnd - lineVectorBegin);
		hullConcave.push_back(temp);
	}

	int posBeginPCL = hullConcave.size();

	for (int i = beginPos; i <= endPos; i++)
	{
		hullConcave.push_back(hullConcaveTemp[i]);
	}

	int posEndPCL = hullConcave.size() - 1;

	lineVectorBegin = hullConcave[hullConcave.size() - 1];
	lineVectorEnd = hullConcave[0];

	for (double i = 0.2; i < 0.9; i += 0.2)
	{
		Point temp = lineVectorBegin + i * (lineVectorEnd - lineVectorBegin);
		hullConcave.push_back(temp);
	}

	/*
		With this the last point ends in an upper corner and it be able to better move the curve top area.
	*/
	std::rotate(hullConcave.begin(), hullConcave.begin() + 1, hullConcave.end());
	int posEndTopArea = hullConcave.size() - 1;
	posBeginTopArea = posBeginTopArea - 1;
	posBeginPCL = posBeginPCL - 1;
	posEndPCL = posEndPCL - 1;

	Point computePoint;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//auto poly = ImplantTools::getPolyLine(hullConcave);
	//std::vector<Point> tExample = { hullConcave[posBeginPCL], hullConcave[posEndPCL], hullConcave[posEndTopArea], hullConcave[posBeginTopArea] };
	//ImplantTools::show(poly, tExample);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (int i = 0; i < hullConcave.size(); i++)
	{
		if (i >= posBeginTopArea && i <= posEndTopArea)
		{
			computePoint = hullConcave[i] + distancePcl * topPlane.getNormalVector();
			finalHull.push_back(computePoint);
		}
		else if (distance != 0)
		{
			if (i + 1 != hullConcave.size())
			{
				computePoint = movePointAtNormal(hullConcave[i], hullConcave[i + 1], myRotation, distance);
				finalHull.push_back(computePoint);
			}
		}
		else
		{
			finalHull.push_back(hullConcave[i]);
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//auto poly2 = ImplantTools::getPolyLine(hullConcave);
	//ImplantTools::show(poly2, finalHull, true);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<Point> concaveSpline = ConvexHull::interpolateSpline(finalHull, amount);
	int initPos = -1;

	for (int i = 0; i < concaveSpline.size(); i++)
	{
		if (obliqueLatPlaneDown.eval(concaveSpline[i]) >= 0 && obliqueMedPlaneDown.eval(concaveSpline[i]) >= 0)
		{
			initPos = i;
			break;
		}
	}

	if (initPos == -1)
	{
		initPos = 0;
	}

	if (initPos > 0)
	{
		std::rotate(concaveSpline.begin(), concaveSpline.begin() + initPos, concaveSpline.end());
	}

	std::vector<Point> finalSplinePointsTKA, finalSplinePointsPKA;

	for (int i = 0; i < concaveSpline.size(); i++)
	{
		if (!(obliqueLatPlaneDown.eval(concaveSpline[i]) >= 0 && obliqueMedPlaneDown.eval(concaveSpline[i]) >= 0))
		{
			finalSplinePointsTKA.push_back(concaveSpline[i]);
		}
	}

	Point initExtreme, endExtreme;

	if (finalSplinePointsTKA.size() > 1)
	{
		initExtreme = finalSplinePointsTKA[0];
		endExtreme = finalSplinePointsTKA[finalSplinePointsTKA.size() - 1];

		cv::Mat initExtremeMat, endExtremeMat;

		initExtremeMat = myRotation * initExtreme.ToMatPoint();
		endExtremeMat = myRotation * endExtreme.ToMatPoint();

		initExtreme = Point(initExtremeMat);
		endExtreme = Point(endExtremeMat);

		if (endExtreme.y < initExtreme.y) 
		{
			std::reverse(finalSplinePointsTKA.begin(), finalSplinePointsTKA.end());
		}
	}
	////////////////////////////////////////////Div PKA section/////////////////////////////////////////////////////////////////////
	//ImplantTools::show(contourMax, finalSplinePointsTKA);
	///////////////////////////////////////////////////////////////////////

	Plane midPlane = myPlane.getPerpendicularPlane(pcl, tubercle);

	if (knee.getSurgerySide() == SurgerySideEnum::KLateral)
	{
		midPlane.movePlane(medPlateau);
		midPlane.reverseByPoint(latPlateau);
		midPlane.movePlane(pcl);
	}
	else
	{
		midPlane.movePlane(latPlateau);
		midPlane.reverseByPoint(medPlateau);
		midPlane.movePlane(pcl);
	}

	Point implantPCL = finalTransformPoint(implant.getPointPCL(), pTransformIn);
	Point implantTuber = finalTransformPoint(implant.getPointTuber(), pTransformIn);
	Point implantSide = finalTransformPoint(implant.getPlateauRefPointDown(), pTransformIn);
	Point implantPlaneSidePoint = finalTransformPoint(implant.getPlaneSidePoint(), pTransformIn);

	Plane midPlaneImplant = myPlane.getPerpendicularPlane(implantPCL, implantTuber);
	midPlaneImplant.reverseByPoint(implantSide);

	Point vectorToPlaneSide = -midPlaneImplant.getNormalVector();
	double sidePlaneThickness = abs(midPlaneImplant.eval(implantPlaneSidePoint)) + distanceSide;

	midPlaneImplant.movePlaneOnNormal(-sidePlaneThickness);

	if (midPlane.eval(finalSplinePointsTKA[0]) < 0)
	{
		for (float i = 0.0; i <= 1.0; i += 0.1)
		{
			Point temp = implantTuber + i * (implantPCL - implantTuber) + sidePlaneThickness * vectorToPlaneSide;
			finalSplinePointsPKA.push_back(myPlane.getProjectionPoint(temp));
		}
	}

	for (int i = 0; i < finalSplinePointsTKA.size(); i++)
	{
		if (midPlaneImplant.eval(finalSplinePointsTKA[i]) > 0)
		{
			finalSplinePointsPKA.push_back(finalSplinePointsTKA[i]);
		}
	}

	if (midPlane.eval(finalSplinePointsTKA[0]) > 0)
	{
		for (float i = 0.5; i <= 1.0; i += 0.1)
		{
			Point temp = implantPCL + i * (implantTuber - implantPCL) + sidePlaneThickness * vectorToPlaneSide;
			finalSplinePointsPKA.push_back(myPlane.getProjectionPoint(temp));
		}
	}

	/////////////////////////////////////////////////////////
	//ImplantTools::show(contourMax, finalSplinePointsPKA);
	//ImplantTools::show(contourMax, finalSplinePointsPKA, true);
	/////////////////////////////////////////////////////////////

	hull = increaseVectorToAmount(finalSplinePointsPKA, amount);

	itk::Vector< double, 3 > translate;
	//ImplantTools::fitEllipse(finalSplinePointsPKA, myPlane.getNormalVector(), hullCenter);
	ConvexHull tempHull(finalSplinePointsPKA, myRotation);
	hullCenter = ImplantTools::getPolygonCenter(tempHull.GetConvexHull(), myPlane.getNormalVector());

	/*auto vectorTest = finalSplinePointsPKA;
	vectorTest.push_back(hullCenter);
	ImplantTools::show(contourMax, vectorTest);*/

	if (tHullSize > 0)
	{
		hullCenter = myPlane.getProjectionPoint(hullCenter);
		Point tTemp = myRotation * (hullCenter.ToMatPoint());
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
