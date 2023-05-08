#include "ImplantTools.hpp"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkKdTreePointLocator.h"
#include "vtkCellLocator.h"
#include <Eigen/Core>
#include <unsupported/Eigen/Splines>
#include <fstream>
#include <opencv2/opencv.hpp>

#include "vtkNamedColors.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSphereSource.h"
#include <itkVersorRigid3DTransform.h>
#include "vtkLine.h"
#include "vtkSelectionNode.h"
#include "vtkSelection.h"
#include "vtkExtractSelectedIds.h"
#include "ConvexHull.hpp"

using namespace THA::IMPLANTS;

void ImplantTools::saveVectorPoints(const std::vector<Point>& pPoints, const std::string& pName, int axis)
{
	ofstream MyFile(pName);
	auto it1 = pPoints.begin();
	auto it2 = pPoints.end();

	for (; it1 != it2; ++it1)
	{
		if (axis == 3)
		{
			MyFile << it1->x << " " << it1->y << " " << it1->z << "\n";
		}
		else if (axis == 2)
		{
			MyFile << it1->x << " " << it1->y << "\n";
		}
		else
		{
			MyFile << it1->x << "\n";
		}
	}

	MyFile.close();
}


cv::Mat ImplantTools::GetRotateZ(const Point& vector)
{
	Point normalXY(0, 0, 1);
	Point rotationAxis = normalXY.cross(vector);
	rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
	double rotationAngle = ImplantTools::getAngleBetweenVectors(normalXY, vector);
	cv::Mat rotation_1 = ImplantTools::getRotateMatrix(rotationAxis, -rotationAngle);
	cv::Mat rotation_2 = ImplantTools::getRotateMatrix(rotationAxis, rotationAngle);
	cv::Mat rotateVector_1 = rotation_1 * vector.ToMatPoint();
	cv::Mat rotateVector_2 = rotation_2 * vector.ToMatPoint();

	cv::Mat rotate;

	double distance_1 = ImplantTools::getAngleBetweenVectors(Point(rotateVector_1), normalXY);
	double distance_2 = ImplantTools::getAngleBetweenVectors(Point(rotateVector_2), normalXY);

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

cv::Mat ImplantTools::GetRotateY(const Point& vector)
{
	Point normalXZ(0, 1, 0);
	Point rotationAxis = normalXZ.cross(vector);
	rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
	double rotationAngle = ImplantTools::getAngleBetweenVectors(normalXZ, vector);
	cv::Mat rotation_1 = ImplantTools::getRotateMatrix(rotationAxis, -rotationAngle);
	cv::Mat rotation_2 = ImplantTools::getRotateMatrix(rotationAxis, rotationAngle);
	cv::Mat rotateVector_1 = rotation_1 * vector.ToMatPoint();
	cv::Mat rotateVector_2 = rotation_2 * vector.ToMatPoint();

	cv::Mat rotate;

	double distance_1 = ImplantTools::getAngleBetweenVectors(Point(rotateVector_1), normalXZ);
	double distance_2 = ImplantTools::getAngleBetweenVectors(Point(rotateVector_2), normalXZ);

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

cv::Mat ImplantTools::GetRotateX(const Point& vector)
{
	Point normalYZ(1, 0, 0);
	Point rotationAxis = normalYZ.cross(vector);
	rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
	double rotationAngle = ImplantTools::getAngleBetweenVectors(normalYZ, vector);
	cv::Mat rotation_1 = ImplantTools::getRotateMatrix(rotationAxis, -rotationAngle);
	cv::Mat rotation_2 = ImplantTools::getRotateMatrix(rotationAxis, rotationAngle);
	cv::Mat rotateVector_1 = rotation_1 * vector.ToMatPoint();
	cv::Mat rotateVector_2 = rotation_2 * vector.ToMatPoint();

	cv::Mat rotate;

	double distance_1 = ImplantTools::getAngleBetweenVectors(Point(rotateVector_1), normalYZ);
	double distance_2 = ImplantTools::getAngleBetweenVectors(Point(rotateVector_2), normalYZ);

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

cv::Mat ImplantTools::getRotateMatrix(const Point& axis, double angle)
{
	cv::Mat rotationMatrix(3, 3, CV_64F);

	Point normaliceAxis = axis;
	normaliceAxis.normalice();

	rotationMatrix.at <double>(0, 0) = 0;
	rotationMatrix.at <double>(1, 0) = normaliceAxis.z;
	rotationMatrix.at <double>(2, 0) = -normaliceAxis.y;

	rotationMatrix.at <double>(0, 1) = -normaliceAxis.z;
	rotationMatrix.at <double>(1, 1) = 0;
	rotationMatrix.at <double>(2, 1) = normaliceAxis.x;

	rotationMatrix.at <double>(0, 2) = normaliceAxis.y;
	rotationMatrix.at <double>(1, 2) = -normaliceAxis.x;
	rotationMatrix.at <double>(2, 2) = 0;

	cv::Mat I = cv::Mat::eye(3, 3, CV_64F);
	cv::Mat result = I + sin(angle)*rotationMatrix + (1.0 - cos(angle))*(rotationMatrix * rotationMatrix);
	return result;
}

cv::Mat ImplantTools::GetGeneralRotateTransformVectors(const Point pFromVector, const Point pToVector)
{
	Point fromVector = pFromVector;
	Point toVector = pToVector;
	fromVector.normalice();
	toVector.normalice();

	Point rotationAxis = fromVector.cross(toVector);
	rotationAxis.normalice();

	double rotationAngle = Line::getAngleBetweenVectors(fromVector, toVector);

	cv::Mat rotation_1 = ImplantTools::getRotateMatrix(rotationAxis, -rotationAngle);
	cv::Mat rotation_2 = ImplantTools::getRotateMatrix(rotationAxis, rotationAngle);

	cv::Mat rotateVector_1 = rotation_1 * fromVector.ToMatPoint();
	cv::Mat rotateVector_2 = rotation_2 * fromVector.ToMatPoint();

	cv::Mat rotate;

	double distance_1 = ImplantTools::getAngleBetweenVectors(Point(rotateVector_1), toVector);
	double distance_2 = ImplantTools::getAngleBetweenVectors(Point(rotateVector_2), toVector);

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

double ImplantTools::getAngleBetweenVectors(const Point& a, const Point& b)
{
	double scalar = a.dot(b);
	double magnitude = sqrt((a.dot(a)) * (b.dot(b)));
	double tCos = scalar / magnitude;
	if (tCos <= -1.0)
	{
		return PI;
	}
	else if (tCos >= 1.0)
	{
		return 0;
	}
	else
	{
		return acos(tCos);
	}
}

double ImplantTools::getAngleBetweenVectorsDegree(const Point& a, const Point& b)
{
	double scalar = a.dot(b);
	double magnitude = sqrt((a.dot(a)) * (b.dot(b)));
	double tCos = scalar / magnitude;
	if (tCos <= -1.0)
	{
		return 180.0;
	}
	else if (tCos >= 1.0)
	{
		return 0;
	}
	else
	{
		return (acos(tCos) * 180.0) / PI;
	}
}

double ImplantTools::getDistanceBetweenPoints(const Point& a, const Point& b, bool square)
{
	Point diff = a - b;
	if (square == false)
	{
		return sqrt(diff.dot(diff));
	}
	else
	{
		return diff.dot(diff);
	}
}

double ImplantTools::getDistanceBetweenPoints(const double a[3], const double b[3], bool square)
{
	Point newA(a[0], a[1], a[2]);
	Point newB(b[0], b[1], b[2]);

	return ImplantTools::getDistanceBetweenPoints(newA, newB, square);
}

Point ImplantTools::getLocalMinimum(const std::vector<Point>& pPoints, const Plane& onPlane, const Point& vectorY, const Plane& nearMinPlane)
{
	cv::Mat rotateZ = ImplantTools::GetRotateZ(onPlane.getNormalVector());
	cv::Mat newVectorY = rotateZ * vectorY.ToMatPoint();
	cv::Mat rotateY = ImplantTools::GetRotateY(Point(newVectorY));

	cv::Mat rotateFinal = rotateY * rotateZ;

	int posMin = -1;
	int posMax = -1;

	int tSize = pPoints.size();
	Point nearMin, pointA, pointB;

	double positiveDistance = 9999999.0;
	double negativeDistance = -9999999.0;
	double tempDistance;

	for (int i = 0; i < tSize; i++)
	{
		tempDistance = nearMinPlane.eval(pPoints[i]);
		if (tempDistance >= 0)
		{
			if (tempDistance < positiveDistance)
			{
				positiveDistance = tempDistance;
				pointA = pPoints[i];
			}
		}
		else
		{
			if (tempDistance > negativeDistance)
			{
				negativeDistance = tempDistance;
				pointB = pPoints[i];
			}
		}
	}

	Line intercept = Line::makeLineWithPoints(pointA, pointB);

	nearMin = nearMinPlane.getInterceptionLinePoint(intercept);

	cv::Mat nearMat = rotateFinal * nearMin.ToMatPoint();
	Point tNear = Point(nearMat);

	double minY = tNear.y;
	double maxY = tNear.y;

	for (int i = 0; i < tSize; i++)
	{
		cv::Mat tempMap = rotateFinal * (pPoints[i].ToMatPoint());
		Point temp = Point(tempMap);

		if (temp.y < minY)
		{
			minY = temp.y;
			posMin = i;
		}

		if (temp.y > maxY)
		{
			maxY = temp.y;
			posMax = i;
		}
	}

	if (abs(tNear.y - minY) < abs(tNear.y - maxY))
	{
		if (posMin != -1)
		{
			return pPoints[posMin];
		}
		else
		{
			return nearMin;
		}

	}
	else
	{
		if (posMax != -1)
		{
			return pPoints[posMax];
		}
		else
		{
			return nearMin;
		}
	}
}

Point ImplantTools::getLocalMinimumTest(const std::list<Point>& pPoints, const Plane& onPlane, const Point& vectorY)
{
	cv::Mat rotateZ = ImplantTools::GetRotateZ(onPlane.getNormalVector());

	cv::Mat newVectorY = rotateZ * vectorY.ToMatPoint();
	cv::Mat rotateY = ImplantTools::GetRotateY(Point(newVectorY));

	cv::Mat rotateFinal = rotateY * rotateZ;

	auto it1 = pPoints.begin();
	auto it2 = pPoints.end();

	double maxX = -9999999;
	double minX = 9999999;

	ofstream MyFile("curve.txt");
	for (; it1 != it2; ++it1)
	{
		cv::Mat pointTempMat = rotateFinal * (*it1).ToMatPoint();
		Point pointTemp = Point(pointTempMat);
		MyFile << pointTemp.x << " " << pointTemp.y << "\n";

		if (pointTemp.x > maxX)
		{
			maxX = pointTemp.x;
		}

		if (pointTemp.x < minX)
		{
			minX = pointTemp.x;
		}
	}

	ofstream MyFileFit("curve_fit.txt");

	//PolyConstraintPoint constraintObj;

	Poly myPoly = polyFit(pPoints, rotateFinal, 7);
	std::cout << "Minx: " << myPoly.getLocalMin() << std::endl;
	for (double i = minX; i <= maxX; i++)
	{
		double y = myPoly.eval(i);

		MyFileFit << i << " " << y << "\n";
	}


	return Point();
}

Point ImplantTools::getLocalMinimum(const std::vector<Point>& pPoints, const Plane& onPlane, const Point& vectorY, int degree)
{
	cv::Mat rotateZ = ImplantTools::GetRotateZ(onPlane.getNormalVector());
	cv::Mat newVectorY = rotateZ * vectorY.ToMatPoint();
	cv::Mat rotateY = ImplantTools::GetRotateY(Point(newVectorY));

	cv::Mat rotateFinal = rotateY * rotateZ;

	PolyConstraintPoint constraintObj;

	Poly myPoly = polyFit(pPoints, rotateFinal, degree, constraintObj);
	cv::Mat resultMat = rotateFinal.inv() * myPoly.getLocalMin().ToMatPoint();
	return Point(resultMat);
}

Point ImplantTools::getLocalMinimum(const std::list<Point>& pPoints, const Plane& onPlane, const Point& vectorY, int degree)
{
	cv::Mat rotateZ = ImplantTools::GetRotateZ(onPlane.getNormalVector());
	cv::Mat newVectorY = rotateZ * vectorY.ToMatPoint();
	cv::Mat rotateY = ImplantTools::GetRotateY(Point(newVectorY));

	cv::Mat rotateFinal = rotateY * rotateZ;

	//PolyConstraintPoint constraintObj;
	Poly myPoly = polyFit(pPoints, rotateFinal, degree);
	cv::Mat resultMat = rotateFinal.inv() * myPoly.getLocalMin().ToMatPoint();
	return Point(resultMat);
}

Point ImplantTools::getLocalMax(const std::list<Point>& pPoints, const Plane& onPlane, const Point& vectorY, int degree)
{
	cv::Mat rotateZ = ImplantTools::GetRotateZ(onPlane.getNormalVector());
	cv::Mat newVectorY = rotateZ * vectorY.ToMatPoint();
	cv::Mat rotateY = ImplantTools::GetRotateY(Point(newVectorY));

	cv::Mat rotateFinal = rotateY * rotateZ;

	//PolyConstraintPoint constraintObj;

	Poly myPoly = polyFit(pPoints, rotateFinal, degree);
	cv::Mat resultMat = rotateFinal.inv() * myPoly.getLocalMax().ToMatPoint();
	return Point(resultMat);
}

Point ImplantTools::getLocalMax(const std::vector<Point>& pPoints, const Plane& onPlane, const Point& vectorY, int degree)
{
	cv::Mat rotateZ = ImplantTools::GetRotateZ(onPlane.getNormalVector());
	cv::Mat newVectorY = rotateZ * vectorY.ToMatPoint();
	cv::Mat rotateY = ImplantTools::GetRotateY(Point(newVectorY));

	cv::Mat rotateFinal = rotateY * rotateZ;

	PolyConstraintPoint constraintObj;

	Poly myPoly = polyFit(pPoints, rotateFinal, degree, constraintObj);
	cv::Mat resultMat = rotateFinal.inv() * myPoly.getLocalMax().ToMatPoint();
	return Point(resultMat);
}

Point ImplantTools::getLocalMinimum(const std::vector<Point>& pPoints, const Plane& onPlane, const Point& vectorY, const Point& nearMin)
{
	cv::Mat rotateZ = ImplantTools::GetRotateZ(onPlane.getNormalVector());
	cv::Mat newVectorY = rotateZ * vectorY.ToMatPoint();
	cv::Mat rotateY = ImplantTools::GetRotateY(Point(newVectorY));

	cv::Mat rotateFinal = rotateY * rotateZ;

	int posMin = -1;
	int posMax = -1;

	cv::Mat nearMat = rotateFinal * nearMin.ToMatPoint();
	Point tNear = Point(nearMat);

	double minY = tNear.y;
	double maxY = tNear.y;

	int tSize = pPoints.size();

	for (int i = 0; i < tSize; i++)
	{
		cv::Mat tempMap = rotateFinal * (pPoints[i].ToMatPoint());
		Point temp = Point(tempMap);

		if (temp.y < minY)
		{
			minY = temp.y;
			posMin = i;
		}

		if (temp.y > maxY)
		{
			maxY = temp.y;
			posMax = i;
		}
	}

	if (abs(tNear.y - minY) < abs(tNear.y - maxY))
	{
		if (posMin != -1)
		{
			return pPoints[posMin];
		}
		else
		{
			return nearMin;
		}

	}
	else
	{
		if (posMax != -1)
		{
			return pPoints[posMax];
		}
		else
		{
			return nearMin;
		}
	}
}

ImplantTools::Poly ImplantTools::polyFit(const std::vector<Point>& pPoints, const cv::Mat& pTransformXY, int order, PolyConstraintPoint constraint)
{
	auto it1 = pPoints.begin();
	auto it2 = pPoints.end();
	std::set<double> filterX;
	std::vector<double> pointsX, pointsY;
	double minX = 9999999999;
	double maxX = -9999999999;

	double Xc = 0;
	double Yc = 0;
	double Z = 0;

	for (; it1 != it2; ++it1)
	{
		cv::Mat pointMat = pTransformXY * (*it1).ToMatPoint();
		Point newPoint = Point(pointMat);
		if (filterX.find(newPoint.x) == filterX.end())
		{
			filterX.insert(newPoint.x);
			Z = newPoint.z;
			pointsX.push_back(newPoint.x);
			pointsY.push_back(newPoint.y);
			Xc += newPoint.x;
			Yc += newPoint.y;

			if (newPoint.x > maxX)
			{
				maxX = newPoint.x;
			}

			if (newPoint.x < minX)
			{
				minX = newPoint.x;
			}
		}
	}

	ImplantTools::Poly obj;
	int tSize = pointsX.size();

	if (tSize < order)
	{
		obj.maxX = 0;
		obj.minX = 0;
		obj.Xc = 0;
		obj.Yc = 0;
		obj.Z = 0;
		obj.isFine = false;
		return obj;
	}

	Xc = Xc / double(tSize);
	Yc = Yc / double(tSize);

	int dependentVarPos = -1;
	double bias;
	std::vector<double> constraintCoeff(order + 1, 0);

	/*
		Substituting two points in the model polynomial and obtaining a new equation to substitute into the system of equations.
	*/

	if (constraint.putConstraid == true)
	{
		cv::Mat pointMatA = pTransformXY * (constraint.a).ToMatPoint();
		cv::Mat pointMatB = pTransformXY * (constraint.b).ToMatPoint();

		Point newPointA = Point(pointMatA);
		Point newPointB = Point(pointMatB);

		bias = newPointA.y - newPointB.y;

		double maxVar = 0;
		for (int i = 1; i <= order; i++)
		{
			constraintCoeff[i] = pow((newPointA.x - Xc), i) - pow((newPointB.x - Xc), i);

			if (abs(constraintCoeff[i]) > abs(maxVar))
			{
				maxVar = constraintCoeff[i];
				dependentVarPos = i;
			}
		}

		if (dependentVarPos > 0)
		{
			for (int i = 1; i <= order; i++)
			{
				constraintCoeff[i] = -(constraintCoeff[i] / maxVar);
			}

			constraintCoeff[dependentVarPos] = 0;

			bias = bias / maxVar;
		}
	}

	int fixOrder = (dependentVarPos > 0) ? order : order + 1;

	Eigen::MatrixXd T(tSize, fixOrder);
	Eigen::VectorXd V(tSize, 1);
	Eigen::VectorXd result;
	int cont = 0;
	double coefVar = 0;

	for (int i = 0; i < tSize; ++i)
	{
		cont = 0;
		if (dependentVarPos > 0)
		{
			coefVar = pow((pointsX.at(i) - Xc), dependentVarPos);
		}

		for (int j = 0; j < order + 1; ++j)
		{
			if (dependentVarPos < 0)
			{
				T(i, j) = pow((pointsX.at(i) - Xc), j);
			}
			else
			{
				if (j != dependentVarPos)
				{
					T(i, cont) = pow((pointsX.at(i) - Xc), j) + (coefVar * constraintCoeff[j]);
					cont++;
				}
			}
		}

		if (dependentVarPos < 0)
		{
			V(i, 0) = (pointsY.at(i) - Yc);
		}
		else
		{
			V(i, 0) = (pointsY.at(i) - Yc) - coefVar * bias;
		}
	}

	result = T.householderQr().solve(V);

	std::vector<double> coeff;
	double dependentVar = bias;

	for (int i = 0; i < fixOrder; i++)
	{
		coeff.push_back(result[i]);

		if (dependentVarPos > 0)
		{
			int pos = 0;
			if (i < dependentVarPos)
			{
				pos = i;
			}
			else
			{
				pos = i + 1;
			}

			dependentVar = dependentVar + constraintCoeff[pos] * coeff[i];
		}
	}

	if (dependentVarPos > 0)
	{
		if (dependentVarPos == coeff.size())
		{
			coeff.push_back(dependentVar);
		}
		else
		{
			coeff.insert(coeff.begin() + dependentVarPos, dependentVar);
		}
	}

	obj.coeff = coeff;
	obj.maxX = maxX;
	obj.minX = minX;
	obj.Xc = Xc;
	obj.Yc = Yc;
	obj.Z = Z;
	obj.isFine = true;

	double diff1 = abs(obj.eval(pointsX[0]) - pointsY[0]);
	double diff2 = abs(obj.eval(pointsX[tSize - 1]) - pointsY[tSize - 1]);
	//std::cout << "Diff: " << diff1 << ";  " << diff2 << " First coef: " << coeff[coeff.size() - 1] << std::endl;

	return obj;
}

ImplantTools::Poly ImplantTools::parabolaFitPCL(const std::vector<Point>& pPoints, const cv::Mat& pTransformXY, const Point& fixPointA, const Point& fixPointB)
{
	ImplantTools::PolyConstraintPoint constraintObj;
	constraintObj.putConstraid = false;

	ImplantTools::Poly tPoly = ImplantTools::polyFit(pPoints, pTransformXY, 6, constraintObj);

	if (tPoly.isFine == false)
	{
		return tPoly;
	}

	Point midPoint = (fixPointA + fixPointB) / 2.;

	cv::Mat pointMat = pTransformXY * (midPoint.ToMatPoint());
	Point temp = Point(pointMat);
	double y = tPoly.eval(temp.x);
	temp.y = y;

	pointMat = pTransformXY.inv() * (temp.ToMatPoint());
	midPoint = Point(pointMat);

	std::vector<Point> parabolaPoints = {fixPointA, fixPointB, midPoint};

	ImplantTools::Poly result = ImplantTools::polyFit(parabolaPoints, pTransformXY, 2, constraintObj);
	return result;
}

ImplantTools::Poly ImplantTools::polyFit(const std::list<Point>& pPoints, const cv::Mat& pTransformXY, int order)
{
	auto it1 = pPoints.begin();
	auto it2 = pPoints.end();
	std::set<double> filterX;
	std::vector<double> pointsX, pointsY;
	double minX = 9999999999;
	double maxX = -9999999999;

	double Xc = 0;
	double Yc = 0;
	double Z = 0;

	for (; it1 != it2; ++it1)
	{
		cv::Mat pointMat = pTransformXY * (*it1).ToMatPoint();
		Point newPoint = Point(pointMat);
		if (filterX.find(newPoint.x) == filterX.end())
		{
			filterX.insert(newPoint.x);
			Z = newPoint.z;
			pointsX.push_back(newPoint.x);
			pointsY.push_back(newPoint.y);
			Xc += newPoint.x;
			Yc += newPoint.y;

			if (newPoint.x > maxX)
			{
				maxX = newPoint.x;
			}

			if (newPoint.x < minX)
			{
				minX = newPoint.x;
			}
		}
	}
	ImplantTools::Poly obj;
	int tSize = pointsX.size();

	if (tSize < order)
	{
		obj.maxX = 0;
		obj.minX = 0;
		obj.Xc = 0;
		obj.Yc = 0;
		obj.Z = 0;
		obj.isFine = false;
		return obj;
	}

	Xc = Xc / double(tSize);
	Yc = Yc / double(tSize);

	Eigen::MatrixXd T(tSize, order + 1);
	Eigen::VectorXd V(tSize, 1);// = Eigen::VectorXd::Map(&pointsY.front(), tSize);
	Eigen::VectorXd result;

	for (int i = 0; i < tSize; ++i)
	{
		for (int j = 0; j < order + 1; ++j)
		{
			T(i, j) = pow((pointsX.at(i) - Xc), j);
		}
		V(i, 0) = (pointsY.at(i) - Yc);
	}

	result = T.householderQr().solve(V);

	std::vector<double> coeff;
	for (int i = 0; i < order + 1; i++)
	{
		coeff.push_back(result[i]);
	}

	obj.coeff = coeff;
	obj.maxX = maxX;
	obj.minX = minX;
	obj.Xc = Xc;
	obj.Yc = Yc;
	obj.Z = Z;
	obj.isFine = true;
	return obj;
}

double ImplantTools::Poly::eval(double x)
{
	if (isFine == false)
	{
		throw std::invalid_argument("The polynomial is undefined.");
	}

	double result = 0;
	x = x - Xc;
	for (int i = 0; i < coeff.size(); i++)
	{
		result += coeff[i] * pow(x, i);
	}
	return (result + Yc);
}

Point ImplantTools::Poly::getLocalMin()
{
	if (isFine == false)
	{
		throw std::invalid_argument("The polynomial is undefined.");
	}

	double step = (abs(maxX - minX)) / 100.0;
	double minValueY = eval(minX);
	double minValueX = minX;
	for (double i = minX; i < maxX; i += step)
	{
		double value = eval(i);
		if (value < minValueY)
		{
			minValueY = value;
			minValueX = i;
		}
	}

	double a = minValueX - step;
	double b = minValueX + step;
	step = (abs(b - a)) / 100.0;

	for (double i = a; i <= b; i += step)
	{
		double value = eval(i);
		if (value < minValueY)
		{
			minValueY = value;
			minValueX = i;
		}
	}

	return Point(minValueX, minValueY, Z);
}

Point ImplantTools::Poly::getLocalMax()
{
	if (isFine == false)
	{
		throw std::invalid_argument("The polynomial is undefined.");
	}

	double step = (abs(maxX - minX)) / 100.0;
	double maxValueY = eval(minX);
	double maxValueX = minX;
	for (double i = minX; i < maxX; i += step)
	{
		double value = eval(i);
		if (value > maxValueY)
		{
			maxValueY = value;
			maxValueX = i;
		}
	}

	double a = maxValueX - step;
	double b = maxValueX + step;
	step = (abs(b - a)) / 100.0;

	for (double i = a; i <= b; i += step)
	{
		double value = eval(i);
		if (value > maxValueY)
		{
			maxValueY = value;
			maxValueX = i;
		}
	}

	return Point(maxValueX, maxValueY, Z);
}

vtkSmartPointer<vtkPolyData> ImplantTools::getMaxContour(const vtkSmartPointer<vtkPolyData> polyData, const Point& pNormal, const Point& pPoint)
{
	vtkNew<vtkPlane> plane;
	vtkNew<vtkCutter> cutter;
	cutter->SetInputData(polyData);

	plane->SetNormal(pNormal.x, pNormal.y, pNormal.z);
	plane->SetOrigin(pPoint.x, pPoint.y, pPoint.z);
	cutter->SetCutFunction(plane);
	cutter->Update();

	auto contour = cutter->GetOutput();

	if (contour->GetNumberOfPoints() > 0)
	{
		vtkNew<vtkPolyDataConnectivityFilter> connectivityFilter;
		connectivityFilter->SetInputData(contour);
		connectivityFilter->SetExtractionModeToLargestRegion();
		connectivityFilter->Update();
		return connectivityFilter->GetOutput();
	}
	else
	{
		return vtkSmartPointer<vtkPolyData>::New();
	}
}

std::vector<std::pair<vtkSmartPointer<vtkPolyData>, vtkSmartPointer<vtkPoints>>> ImplantTools::getAllContours(const vtkSmartPointer<vtkPolyData> polyData, const Point& pNormal, const Point& pPoint)
{
	std::vector<std::pair<vtkSmartPointer<vtkPolyData>, vtkSmartPointer<vtkPoints>>> result;
	vtkNew<vtkPlane> plane;
	vtkNew<vtkCutter> cutter;
	cutter->SetInputData(polyData);

	plane->SetNormal(pNormal.x, pNormal.y, pNormal.z);
	plane->SetOrigin(pPoint.x, pPoint.y, pPoint.z);
	cutter->SetCutFunction(plane);
	cutter->Update();

	auto contour = cutter->GetOutput();

	if (contour->GetNumberOfPoints() > 0)
	{
		vtkNew<vtkPolyDataConnectivityFilter> connectivityFilter;
		connectivityFilter->SetInputData(contour);

		connectivityFilter->SetExtractionModeToAllRegions();
		connectivityFilter->Update();
		int connected_components = connectivityFilter->GetNumberOfExtractedRegions();

		for (int i = 0; i < connected_components; i++)
		{
			vtkSmartPointer<vtkPolyData> region = vtkSmartPointer<vtkPolyData>::New();
			connectivityFilter->SetExtractionModeToSpecifiedRegions();
			connectivityFilter->AddSpecifiedRegion(i);
			connectivityFilter->Update();
			region->ShallowCopy(connectivityFilter->GetOutput());

			vtkSmartPointer<vtkIdTypeArray> regionPointIds = vtkSmartPointer<vtkIdTypeArray>::New();
			regionPointIds->SetNumberOfComponents(1);
			for (vtkIdType j = 0; j < region->GetNumberOfCells(); j++)
			{
				vtkCell* cell = region->GetCell(j);
				for (vtkIdType k = 0; k < cell->GetNumberOfPoints(); k++)
				{
					vtkIdType pointId = cell->GetPointId(k);
					regionPointIds->InsertNextValue(pointId);
				}
			}

			// Seleccionar solo los puntos que pertenecen a la región actual
			vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
			selectionNode->SetFieldType(vtkSelectionNode::POINT);
			selectionNode->SetContentType(vtkSelectionNode::INDICES);
			selectionNode->SetSelectionList(regionPointIds);

			vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
			selection->AddNode(selectionNode);

			vtkSmartPointer<vtkExtractSelectedIds> extractSelectedIds = vtkSmartPointer<vtkExtractSelectedIds>::New();
			extractSelectedIds->SetInputData(0, region);
			extractSelectedIds->SetInputData(1, selection);
			extractSelectedIds->Update();

			vtkSmartPointer<vtkPolyData> regionPoints = vtkSmartPointer<vtkPolyData>::New();
			regionPoints->ShallowCopy(extractSelectedIds->GetOutput());

			result.push_back(std::make_pair(region, regionPoints->GetPoints()));

			connectivityFilter->DeleteSpecifiedRegion(i);
		}
	}
	return result;
}

vtkSmartPointer<vtkPolyData> ImplantTools::getContours(const vtkSmartPointer<vtkPolyData> polyData, const Point& pNormal, const Point& pPoint)
{
	vtkNew<vtkPlane> plane;
	vtkNew<vtkCutter> cutter;
	cutter->SetInputData(polyData);

	plane->SetNormal(pNormal.x, pNormal.y, pNormal.z);
	plane->SetOrigin(pPoint.x, pPoint.y, pPoint.z);
	cutter->SetCutFunction(plane);
	cutter->Update();

	auto contour = cutter->GetOutput();

	return contour;
}

void ImplantTools::ExtractSortLines(const vtkSmartPointer<vtkPolyData> polyData, std::list<std::pair<vtkIdType, vtkIdType>>& lines)
{
	vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
	if (points->GetNumberOfPoints() == 0)
	{
		return;
	}
	vtkSmartPointer<vtkCellArray> cells = polyData->GetLines();
	cells->InitTraversal();
	vtkNew<vtkIdList> myList;
	std::list<std::pair<vtkIdType, vtkIdType>> restLines;

	while (cells->GetNextCell(myList))
	{
		vtkIdType id1 = myList->GetId(0);
		vtkIdType id2 = myList->GetId(1);

		if (lines.size() == 0)
		{
			lines.push_back(std::make_pair(id1, id2));
		}
		else
		{
			std::pair<vtkIdType, vtkIdType> front = lines.front();
			std::pair<vtkIdType, vtkIdType> back = lines.back();

			if (id1 == front.first)
			{
				lines.push_front(std::make_pair(id2, id1));
			}
			else if (id2 == front.first)
			{
				lines.push_front(std::make_pair(id1, id2));
			}
			else if (id1 == back.second)
			{
				lines.push_back(std::make_pair(id1, id2));
			}
			else if (id2 == back.second)
			{
				lines.push_back(std::make_pair(id2, id1));
			}
			else
			{
				restLines.push_back(std::make_pair(id1, id2));
			}
		}

	}

	vtkIdType size1 = restLines.size();

	for (int i = 0; i < size1; i++)
	{
		auto it1 = restLines.begin();
		auto it2 = restLines.end();
		for (; it1 != it2; ++it1)
		{
			vtkIdType id1 = (*it1).first;
			vtkIdType id2 = (*it1).second;

			std::pair<vtkIdType, vtkIdType> front = lines.front();
			std::pair<vtkIdType, vtkIdType> back = lines.back();

			if (id1 == front.first)
			{
				lines.push_front(std::make_pair(id2, id1));
				restLines.erase(it1);
				break;
			}
			else if (id2 == front.first)
			{
				lines.push_front(std::make_pair(id1, id2));
				restLines.erase(it1);
				break;
			}
			else if (id1 == back.second)
			{
				lines.push_back(std::make_pair(id1, id2));
				restLines.erase(it1);
				break;
			}
			else if (id2 == back.second)
			{
				lines.push_back(std::make_pair(id2, id1));
				restLines.erase(it1);
				break;
			}
		}
	}
}

double ImplantTools::GetInterceptionWithLine(const vtkSmartPointer<vtkImplicitPolyDataDistance> polyDistance, const Point& pRefPoint, const Point& p2, Point& result)
{
	Point vector = pRefPoint - p2;
	vector.normalice();

	double error = 0.1;
	int iter = 20;
	double currentDistance = 100000;
	result = pRefPoint;
	double myClosest[3];

	do
	{
		double pnt[3] = { result.x, result.y, result.z };
		currentDistance = abs(polyDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest));

		Point a = result + currentDistance * vector;
		Point b = result - currentDistance * vector;

		double pntA[3] = { a.x, a.y, a.z };
		double pntB[3] = { b.x, b.y, b.z };

		double distanceA = abs(polyDistance->EvaluateFunctionAndGetClosestPoint(pntA, myClosest));
		double distanceB = abs(polyDistance->EvaluateFunctionAndGetClosestPoint(pntB, myClosest));

		if (distanceA < distanceB)
		{
			result = Point(pntA[0], pntA[1], pntA[2]);
			currentDistance = distanceA;
		}
		else
		{
			result = Point(pntB[0], pntB[1], pntB[2]);
			currentDistance = distanceB;
		}
		iter--;

	} while (currentDistance > error && iter > 0);

	return currentDistance;
}

bool ImplantTools::GetInterceptionWithSegment(const vtkSmartPointer<vtkPolyData> poly, const Point& p1, const Point& p2, Point& result)
{
	vtkNew<vtkCellLocator> cellLocator;
	cellLocator->SetDataSet(poly);
	cellLocator->BuildLocator();

	double a[3] = { p1.x, p1.y, p1.z };
	double b[3] = { p2.x, p2.y, p2.z };

	double tol = 0;
	double t;
	double x[3];
	double pcoords[3];
	int subId;

	int intercept = cellLocator->IntersectWithLine(a, b, tol, t, x, pcoords, subId);
	if (intercept == 0)
	{
		return false;
	}
	else
	{
		result.x = x[0];
		result.y = x[1];
		result.z = x[2];
		return true;
	}
}

vtkIdType ImplantTools::GetNearestPoints(const vtkSmartPointer<vtkPolyData> poly, const Point& pPoint)
{
	vtkNew<vtkKdTreePointLocator> kDTree;
	kDTree->SetDataSet(poly);
	kDTree->BuildLocator();

	double pnt[3] = { pPoint.x, pPoint.y, pPoint.z };

	vtkIdType nearPoint = kDTree->FindClosestPoint(pnt);

	return nearPoint;
}

std::pair<double, Point> ImplantTools::GetDistancePlaneToSurface(const vtkSmartPointer<vtkPolyData> poly, const Plane& pPlane, const Plane& pUseOneSide)
{
	double distance = -1;
	double temp;
	Point closest;
	vtkSmartPointer<vtkPoints> points = poly->GetPoints();
	vtkIdType tSize = points->GetNumberOfPoints();

	for (int i = 0; i < tSize; i++)
	{
		double pnt[3];
		points->GetPoint(i, pnt);

		if (pUseOneSide.getIsInit() == true)
		{
			if (pUseOneSide.eval(pnt) < 0)
			{
				continue;
			}
		}

		temp = pPlane.getDistanceFromPoint(Point(pnt[0], pnt[1], pnt[2]));

		if (distance < 0 || distance > temp)
		{
			distance = temp;
			closest = Point(pnt[0], pnt[1], pnt[2]);
		}
	}

	return std::make_pair(distance, closest);
}

std::pair<int, float> ImplantTools::GetNearestPointToLine(const std::vector<Point>& pPoints, const Line& pLine, const Plane& pConditionPlane)
{
	int tSize = pPoints.size();
	int pos = -1;
	float distance = -1;
	for (int i = 0; i < tSize; i++)
	{
		if (pConditionPlane.getIsInit() == true)
		{
			if (pConditionPlane.eval(pPoints[i]) < 0)
			{
				continue;
			}
		}

		float temp = pLine.getDistanceFromPoint(pPoints[i]);
		if (temp < distance || distance < 0)
		{
			distance = temp;
			pos = i;
		}
	}

	return std::make_pair(pos, distance);
}

vtkIdType ImplantTools::GetNearestPoints(const vtkSmartPointer<vtkPolyData> poly, const Line& pLine, const Plane& pPlane1, const Plane& pPlane2)
{
	double distance = -1;
	double temp;
	vtkSmartPointer<vtkPoints> points = poly->GetPoints();
	vtkIdType tSize = points->GetNumberOfPoints();
	vtkIdType result = 0;
	for (int i = 0; i < tSize; i++)
	{
		double pnt[3];
		points->GetPoint(i, pnt);

		if (pPlane1.getIsInit() == true)
		{
			if (pPlane1.eval(pnt) < 0)
			{
				continue;
			}
		}

		if (pPlane2.getIsInit() == true)
		{
			if (pPlane2.eval(pnt) < 0)
			{
				continue;
			}
		}

		temp = pLine.getDistanceFromPoint(Point(pnt[0], pnt[1], pnt[2]));

		if (distance < 0 || distance > temp)
		{
			distance = temp;
			result = i;
		}
	}
	return result;
}

Point ImplantTools::GetFarestPoint(const vtkSmartPointer<vtkPolyData> poly, const Plane& pPlane)
{
	vtkSmartPointer<vtkPoints> points = poly->GetPoints();
	vtkIdType tSize = points->GetNumberOfPoints();
	double distance = 0;
	Point result;
	for (vtkIdType i = 0; i < tSize; i++)
	{
		double pnt[3];
		points->GetPoint(i, pnt);
		double temp = pPlane.eval(Point(pnt[0], pnt[1], pnt[2]));
		if (i == 0 || temp > distance)
		{
			distance = temp;
			result = Point(pnt[0], pnt[1], pnt[2]);
		}
	}
	return result;
}

void ImplantTools::GetPointsOnContourSort(const vtkSmartPointer<vtkPolyData> contour, const std::vector<Plane>& pCondition, const Point& nearPoint, std::vector<Point>& pResult)
{
	std::vector<vtkIdType> allPoints;

	std::list<std::pair<vtkIdType, vtkIdType>> lines;
	ImplantTools::ExtractSortLines(contour, lines);

	std::list<std::pair<vtkIdType, vtkIdType>>::iterator it11, it22;
	it11 = lines.begin();
	it22 = lines.end();

	double distance = -1;
	int pos = 0;
	bool isFine;

	for (; it11 != it22; ++it11)
	{
		double pnt[3];
		contour->GetPoint(it11->first, pnt);

		isFine = true;
		for (int j = 0; j < pCondition.size(); j++)
		{
			if (pCondition[j].eval(pnt) <= 0)
			{
				isFine = false;
				break;
			}
		}

		if (isFine == true)
		{
			Point tempPoint = Point(pnt[0], pnt[1], pnt[2]);
			double temp = getDistanceBetweenPoints(tempPoint, nearPoint);
			if (distance < 0 || distance > temp)
			{
				distance = temp;
				pos = allPoints.size();
			}
		}

		allPoints.push_back(it11->first);
	}

	std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());

	isFine = true;
	for (int j = 0; j < pCondition.size(); j++)
	{
		double pnt[3];
		contour->GetPoint(allPoints[1], pnt);

		if (pCondition[j].eval(pnt) <= 0)
		{
			isFine = false;
			break;
		}
	}

	if (isFine == false)
	{
		vtkIdType myTemp = allPoints[0];
		allPoints.push_back(myTemp);
		allPoints[0] = allPoints[1];
		std::reverse(allPoints.begin(), allPoints.end());
	}

	vtkIdType tSize = allPoints.size();


	for (vtkIdType i = 0; i < tSize; i++)
	{
		double pnt[3];
		contour->GetPoint(allPoints[i], pnt);
		isFine = true;
		for (int j = 0; j < pCondition.size(); j++)
		{
			if (pCondition[j].eval(pnt) <= 0)
			{
				isFine = false;
				break;
			}
		}

		if (isFine == true)
		{
			pResult.push_back(Point(pnt[0], pnt[1], pnt[2]));
		}
	}
}

void ImplantTools::GetPointsOnContourSort(const vtkSmartPointer<vtkPolyData> contour, const std::vector<Plane>& pCondition, const Plane& nearPlane, std::vector<Point>& pResult)
{
	std::vector<vtkIdType> allPoints;

	std::list<std::pair<vtkIdType, vtkIdType>> lines;
	ImplantTools::ExtractSortLines(contour, lines);

	std::list<std::pair<vtkIdType, vtkIdType>>::iterator it11, it22;
	it11 = lines.begin();
	it22 = lines.end();

	double distance = -1;
	int pos = 0;
	bool isFine;

	for (; it11 != it22; ++it11)
	{
		double pnt[3];
		contour->GetPoint(it11->first, pnt);

		isFine = true;
		for (int j = 0; j < pCondition.size(); j++)
		{
			if (pCondition[j].eval(pnt) <= 0)
			{
				isFine = false;
				break;
			}
		}

		if (isFine == true)
		{
			double temp = abs(nearPlane.eval(pnt));
			if (distance < 0 || distance > temp)
			{
				distance = temp;
				pos = allPoints.size();
			}
		}

		allPoints.push_back(it11->first);
	}

	std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());

	isFine = true;
	for (int j = 0; j < pCondition.size(); j++)
	{
		double pnt[3];
		contour->GetPoint(allPoints[1], pnt);

		if (pCondition[j].eval(pnt) <= 0)
		{
			isFine = false;
			break;
		}
	}

	if (isFine == false)
	{
		vtkIdType myTemp = allPoints[0];
		allPoints.push_back(myTemp);
		allPoints[0] = allPoints[1];
		std::reverse(allPoints.begin(), allPoints.end());
	}

	vtkIdType tSize = allPoints.size();

	for (vtkIdType i = 0; i < tSize; i++)
	{
		double pnt[3];
		contour->GetPoint(allPoints[i], pnt);

		isFine = true;
		for (int j = 0; j < pCondition.size(); j++)
		{
			if (pCondition[j].eval(pnt) <= 0)
			{
				isFine = false;
				break;
			}
		}

		if (isFine == true)
		{
			pResult.push_back(Point(pnt[0], pnt[1], pnt[2]));
		}
	}
}

Point ImplantTools::GetFarestPoint(const vtkSmartPointer<vtkPolyData> poly, const Plane& pPlane, const Plane& pCondition)
{
	vtkSmartPointer<vtkPoints> points = poly->GetPoints();
	vtkIdType tSize = points->GetNumberOfPoints();
	double distance = 0;
	Point result;
	for (vtkIdType i = 0; i < tSize; i++)
	{
		double pnt[3];
		points->GetPoint(i, pnt);

		if (pCondition.eval(pnt) < 0)
		{
			continue;
		}

		double temp = pPlane.eval(Point(pnt[0], pnt[1], pnt[2]));
		if (i == 0 || temp > distance)
		{
			distance = temp;
			result = Point(pnt[0], pnt[1], pnt[2]);
		}
	}
	return result;
}

bool ImplantTools::isPointInsideSphere(const Point& center, double radius, const Point& pPoint)
{
	double dist = ImplantTools::getDistanceBetweenPoints(center, pPoint);
	if (dist <= radius)
	{
		return true;
	}
	else
	{
		return false;
	}
}

std::pair<cv::Point2d, double> ImplantTools::findCircle(const cv::Point2d& p1, const cv::Point2d& p2, const cv::Point2d& p3)
{
	double x12 = p1.x - p2.x;
	double x13 = p1.x - p3.x;

	double y12 = p1.y - p2.y;
	double y13 = p1.y - p3.y;

	double y31 = p3.y - p1.y;
	double y21 = p2.y - p1.y;

	double x31 = p3.x - p1.x;
	double x21 = p2.x - p1.x;

	// x1^2 - x3^2
	double sx13 = pow(p1.x, 2) - pow(p3.x, 2);

	// y1^2 - y3^2
	double sy13 = pow(p1.y, 2) - pow(p3.y, 2);

	double sx21 = pow(p2.x, 2) - pow(p1.x, 2);
	double sy21 = pow(p2.y, 2) - pow(p1.y, 2);

	double f = ((sx13) * (x12)
		+(sy13) * (x12)
		+(sx21) * (x13)
		+(sy21) * (x13))
		/ (2 * ((y31) * (x12)-(y21) * (x13)));
	double g = ((sx13) * (y12)
		+(sy13) * (y12)
		+(sx21) * (y13)
		+(sy21) * (y13))
		/ (2 * ((x31) * (y12)-(x21) * (y13)));

	int c = -pow(p1.x, 2) - pow(p1.y, 2) - 2 * g * p1.x - 2 * f * p1.y;

	// eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0
	// where centre is (h = -g, k = -f) and radius r
	// as r^2 = h^2 + k^2 - c
	double h = -g;
	double k = -f;
	double sqr_of_r = h * h + k * k - c;

	// r is the radius
	double r = sqrt(sqr_of_r);

	return std::make_pair(cv::Point2d(h, k), r);
}

std::pair<Point, double> ImplantTools::fitSphere(const std::vector<Point>& pPoints)
{
	if (pPoints.size() < 4)
	{
		return std::make_pair(Point(), -1);
	}

	std::vector<Point>::const_iterator it1, it2;
	it1 = pPoints.begin();
	it2 = pPoints.end();
	std::vector<cv::Mat> points;
	std::vector<double> bias;

	cv::Mat A(pPoints.size(), 4, CV_64F);
	int cont = 0;
	for (; it1 != it2; ++it1)
	{
		Point temp = (*it1);
		A.at<double>(cont, 0) = temp.x;
		A.at<double>(cont, 1) = temp.y;
		A.at<double>(cont, 2) = temp.z;
		A.at<double>(cont, 3) = 1.0;
		bias.push_back(temp.dot(temp));
		cont++;
	}

	cv::Mat B = cv::Mat(bias.size(), 1, CV_64F, bias.data());

	cv::Mat dx;

	bool result = cv::solve(A, B, dx, cv::DECOMP_SVD);

	if (result == false)
	{
		return std::make_pair(Point(), -1);
	}

	Point center = (Point(dx.at<double>(0, 0), dx.at<double>(1, 0), dx.at<double>(2, 0))) / 2.0;
	double radius = sqrt(center.dot(center) + dx.at<double>(3, 0));
	return std::make_pair(center, radius);
}

std::pair<double, double> ImplantTools::fitEllipse(const vtkSmartPointer<vtkPolyData> pContour, const Point& pNormal, Point& center)
{
	cv::Mat rotate = GetRotateZ(pNormal);
	vtkSmartPointer<vtkPoints> pointList = pContour->GetPoints();
	int tSize = pointList->GetNumberOfPoints();

	std::vector<cv::Point2f> coplanar2d;
	double z = 0;

	for (int i = 0; i < tSize; i++)
	{
		double pnt[3];
		pointList->GetPoint(i, pnt);
		Point newPoint(pnt[0], pnt[1], pnt[2]);
		cv::Mat transformPointMat = rotate * newPoint.ToMatPoint();
		Point tPoint = Point(transformPointMat);
		coplanar2d.push_back(cv::Point2f(tPoint.x, tPoint.y));
		z = tPoint.z;
	}

	cv::RotatedRect box = cv::fitEllipse(coplanar2d);
	std::pair<double, double> result = std::make_pair(box.size.width, box.size.height);

	Point newCenter = Point(box.center.x, box.center.y, z);
	cv::Mat resultCenter = rotate.inv() * newCenter.ToMatPoint();
	center = Point(resultCenter);

	return result;
}

int ImplantTools::orientation(const cv::Point2d& p1, const cv::Point2d& p2, const cv::Point2d& p3)
{
	double val = (p2.y - p1.y) * (p3.x - p2.x) - (p2.x - p1.x) * (p3.y - p2.y);

	if (val == 0) return 0;  // collinear

	return (val > 0) ? 1 : 2; // clock or counterclock wise
}

int ImplantTools::orientation(const Point& p1, const Point& p2, const Point& p3, const cv::Mat& pRotationZ)
{
	cv::Mat aMat = pRotationZ * p1.ToMatPoint();
	cv::Mat bMat = pRotationZ * p2.ToMatPoint();
	cv::Mat cMat = pRotationZ * p3.ToMatPoint();

	Point a = Point(aMat);
	Point b = Point(bMat);
	Point c = Point(cMat);

	return orientation(cv::Point2d(a.x, a.y), cv::Point2d(b.x, b.y), cv::Point2d(c.x, c.y));
}

int ImplantTools::GetCornerPointOnContour(const std::vector<Point>& pPoints, const Point& pCenter, const Point& pExtremeA, const Point& pExtremeB, Plane pConditionPlane, double pWeightA, double pWeightB)
{
	Point vectorA = pExtremeA - pCenter;
	Point vectorB = pExtremeB - pCenter;
	double angle = getAngleBetweenVectorsDegree(vectorA, vectorB);

	if (angle <= 1 || angle >= 179)
	{
		return -1;
	}

	Plane planeA, planeB;
	planeA.init(vectorB, pCenter);
	planeB.init(vectorA, pCenter);

	planeA.reverseByPoint(pExtremeB);
	planeB.reverseByPoint(pExtremeA);

	double tResize = 1000 * vectorA.dot(vectorA) + 1000 * vectorB.dot(vectorB);

	vectorA.normalice();
	vectorB.normalice();

	Point farA = pCenter + pWeightA * tResize * vectorA;
	Point farB = pCenter + pWeightB * tResize * vectorB;

	Line myLine = Line::makeLineWithPoints(farA, farB);

	int tSize = pPoints.size();
	double distance = -1;
	int pos = -1;

	for (int i = 0; i < tSize; i++)
	{
		if (pConditionPlane.getIsInit() == true)
		{
			if (pConditionPlane.eval(pPoints[i]) < 0)
			{
				continue;
			}
		}

		if (planeA.eval(pPoints[i]) > 0 && planeB.eval(pPoints[i]) > 0)
		{
			double tDist = myLine.getDistanceFromPoint(pPoints[i]);

			if (distance < 0 || tDist < distance)
			{
				distance = tDist;
				pos = i;
			}
		}
	}
	return pos;
}

cv::Mat ImplantTools::Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer transform)
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

cv::Mat ImplantTools::Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer transform)
{
	itk::Vector< double, 3 > translate = transform->GetOffset();
	Point temp(translate[0], translate[1], translate[2]);
	return temp.ToMatPoint();
}

cv::Mat ImplantTools::Rigid3DTransformToCV(const itk::Rigid3DTransform<>::Pointer transform)
{
	itk::Matrix< double, 3, 3 > rotation = transform->GetMatrix();
	itk::Vector< double, 3 > translate = transform->GetOffset();
	double* matrix = new double[16];

	matrix[0] = rotation[0][0];
	matrix[1] = rotation[0][1];
	matrix[2] = rotation[0][2];
	matrix[3] = translate[0];

	matrix[4] = rotation[1][0];
	matrix[5] = rotation[1][1];
	matrix[6] = rotation[1][2];
	matrix[7] = translate[1];

	matrix[8] = rotation[2][0];
	matrix[9] = rotation[2][1];
	matrix[10] = rotation[2][2];
	matrix[11] = translate[2];

	matrix[12] = 0;
	matrix[13] = 0;
	matrix[14] = 0;
	matrix[15] = 1;

	cv::Mat matrixCV(4, 4, CV_64FC1, matrix);
	return matrixCV;
}

cv::Mat ImplantTools::JoinRigidTransform(const cv::Mat& rotation, const cv::Mat& translation)
{
	double* matrix = new double[16];

	matrix[0] = rotation.at<double>(0, 0);
	matrix[1] = rotation.at<double>(0, 1);
	matrix[2] = rotation.at<double>(0, 2);
	matrix[3] = translation.at<double>(0, 0);

	matrix[4] = rotation.at<double>(1, 0);
	matrix[5] = rotation.at<double>(1, 1);
	matrix[6] = rotation.at<double>(1, 2);
	matrix[7] = translation.at<double>(1, 0);

	matrix[8] = rotation.at<double>(2, 0);
	matrix[9] = rotation.at<double>(2, 1);
	matrix[10] = rotation.at<double>(2, 2);
	matrix[11] = translation.at<double>(2, 0);

	matrix[12] = 0;
	matrix[13] = 0;
	matrix[14] = 0;
	matrix[15] = 1;

	cv::Mat matrixCV(4, 4, CV_64FC1, matrix);
	return matrixCV;
}

itk::Rigid3DTransform<>::Pointer ImplantTools::getITKTransformFromCV(const cv::Mat& fullTransform)
{
	itk::Matrix< double, 3, 3 > rotation;
	itk::Vector< double, 3 > translate;

	rotation[0][0] = fullTransform.at<double>(0, 0);
	rotation[0][1] = fullTransform.at<double>(0, 1);
	rotation[0][2] = fullTransform.at<double>(0, 2);

	rotation[1][0] = fullTransform.at<double>(1, 0);
	rotation[1][1] = fullTransform.at<double>(1, 1);
	rotation[1][2] = fullTransform.at<double>(1, 2);

	rotation[2][0] = fullTransform.at<double>(2, 0);
	rotation[2][1] = fullTransform.at<double>(2, 1);
	rotation[2][2] = fullTransform.at<double>(2, 2);

	translate[0] = fullTransform.at<double>(0, 3);
	translate[1] = fullTransform.at<double>(1, 3);
	translate[2] = fullTransform.at<double>(2, 3);

	itk::Rigid3DTransform<double>::Pointer transform = itk::VersorRigid3DTransform<double>::New();

	transform->SetMatrix(rotation);
	transform->SetOffset(translate);

	return transform;
}

void ImplantTools::replaceVectorRange(const std::vector<Point>& pNewValues, int pInitPos, int pEndPos, std::vector<Point>& pVectorToReplace, bool sameOrder, bool includeInitPos)
{

	if (pEndPos < pInitPos || pVectorToReplace.size() <= pEndPos || pInitPos < 0)
	{
		return;
	}

	std::vector<Point> temp = pNewValues;
	if (sameOrder == false && temp.size() > 1)
	{
		std::reverse(temp.begin(), temp.end());
	}

	std::vector<Point> tempJoin;

	int tSizeInit = pInitPos;
	if (includeInitPos == false)
	{
		tSizeInit = tSizeInit - 1;
	}

	for (int i = 0; i <= tSizeInit; i++)
	{
		tempJoin.push_back(pVectorToReplace[i]);
	}

	for (int i = 0; i < temp.size(); i++)
	{
		tempJoin.push_back(temp[i]);
	}

	if (pInitPos == pEndPos)
	{
		pEndPos = pEndPos + 1;
	}

	for (int i = pEndPos; i < pVectorToReplace.size(); i++)
	{
		tempJoin.push_back(pVectorToReplace[i]);
	}

	pVectorToReplace.clear();
	pVectorToReplace = tempJoin;
}

void ImplantTools::squareCorner(const Point& pBasePoint, const Point& pCenter, const Point& pSecondPoint, const Plane& pRightSide, std::vector<Point>& pData, double pSlopeAngle)
{
	Plane basePlane;
	Point vector1 = pBasePoint - pCenter;
	Point vector2 = pSecondPoint - pCenter;
	vector1.normalice();
	vector2.normalice();

	Point farBase = pBasePoint + 1000 * vector1;
	Point farSecond = pSecondPoint + 1000 * vector2;

	basePlane.init(vector1, farBase);
	Line myLine(vector1, farSecond);

	////////////////////////////////////////////Find extreme points

	auto it1 = pData.begin();
	auto it2 = pData.end();

	double distance, baseDist = -1, secondDist = -1;

	Point tBasePoint, tSecondPoint;

	int basePos, secondPos, cont = 0;

	for (; it1 != it2; it1++)
	{
		if (pRightSide.eval(*it1) < 0)
		{
			cont++;
			continue;
		}

		distance = basePlane.getDistanceFromPoint(*it1);

		if (distance < baseDist || baseDist == -1)
		{
			baseDist = distance;
			tBasePoint = *it1;
			basePos = cont;
		}

		distance = myLine.getDistanceFromPoint(*it1);

		if (distance < secondDist || secondDist == -1)
		{
			secondDist = distance;
			tSecondPoint = *it1;
			secondPos = cont;
		}

		cont++;
	}

	bool tIsAreaInsidePos;

	//////////////////////////////////////////
	/*
	looking for the area with the greatest number of points to avoid false positives with the areas of belonging between the planes
	*/

	int distInside, distBegin, distEnd;
	distInside = abs(basePos - secondPos);

	if (basePos < secondPos)
	{
		distBegin = basePos;
		distEnd = pData.size() - secondPos;
	}

	else
	{
		distBegin = secondPos;
		distEnd = pData.size() - basePos;
	}

	//////////////////////////////////// Detecting where the area to square is in the pData vector.

	Plane currentPlane, pBase, pSecond;
	currentPlane.init(pBasePoint, pCenter, pSecondPoint);
	pBase = currentPlane.getPerpendicularPlane(tSecondPoint, pCenter);
	pSecond = currentPlane.getPerpendicularPlane(tBasePoint, pCenter);

	pBase.reverseByPoint(farBase);
	pSecond.reverseByPoint(farSecond);

	if (distInside >= distBegin && distBegin >= distEnd)
	{
		int evalPos = (basePos + secondPos) / 2;

		if (pBase.eval(pData[evalPos]) > 0 && pSecond.eval(pData[evalPos]) > 0)
		{
			tIsAreaInsidePos = true;
		}
		else
		{
			tIsAreaInsidePos = false;
		}
	}
	else
	{
		int evalPos;
		int mid = (distEnd + distBegin) / 2;

		if (mid < distEnd)
		{
			evalPos = pData.size() - distEnd + mid;
		}
		else
		{
			evalPos = mid - distEnd;
		}

		if (pBase.eval(pData[evalPos]) > 0 && pSecond.eval(pData[evalPos]) > 0)
		{
			tIsAreaInsidePos = false;
		}
		else
		{
			tIsAreaInsidePos = true;
		}
	}

	/////////////////////////////////////////// Generate the new corner
	basePlane.movePlane(tBasePoint);
	myLine.setPoint(tSecondPoint);
	Point intercepTop = basePlane.getInterceptionLinePoint(myLine);
	Point lastTopPoint = tBasePoint;
	vector1 = tSecondPoint - intercepTop;
	vector2 = lastTopPoint - intercepTop;
	double angle = ImplantTools::getAngleBetweenVectorsDegree(vector1, vector2);

	Point newPoint = tSecondPoint;
	int posNewPoint = secondPos;

	if (angle < pSlopeAngle - 0.5)
	{
		//applying the law of sines.
		double myAngle = 180. - pSlopeAngle;
		double dist1 = sqrt(vector1.dot(vector1));
		double dist2 = sqrt(vector2.dot(vector2));
		double angle2 = 180 - (angle + myAngle);
		double moveDist = (dist1 / sin(myAngle * PI / 180.)) * sin(angle2 * PI / 180.);

		vector2.normalice();
		if (moveDist < dist2)
		{
			newPoint = intercepTop + moveDist * vector2;
		}
		else
		{
			newPoint = lastTopPoint;
		}

		myLine.setDirectVector((newPoint - tSecondPoint));
		myLine.setPoint(intercepTop);

		it1 = pData.begin();
		it2 = pData.end();
		secondDist = -1;

		cont = 0;
		for (; it1 != it2; ++it1)
		{
			if (pRightSide.eval(*it1) < 0)
			{
				cont++;
				continue;
			}

			distance = myLine.getDistanceFromPoint(*it1);
			if (distance < secondDist || secondDist < 0)
			{
				secondDist = distance;
				newPoint = *it1;
				posNewPoint = cont;
			}
			cont++;
		}

		myLine.setPoint(newPoint);
		intercepTop = basePlane.getInterceptionLinePoint(myLine);
	}

	std::vector<Point> squarePoints, tempAll;

	if (intercepTop != lastTopPoint)
	{
		for (double i = 0.2; i < 0.9; i += 0.2)
		{
			Point temp = lastTopPoint + i * (intercepTop - lastTopPoint);
			squarePoints.push_back(temp);
		}
		squarePoints.push_back(intercepTop);
	}

	for (double i = 0.2; i < 0.9; i += 0.2)
	{
		Point temp = intercepTop + i * (newPoint - intercepTop);
		squarePoints.push_back(temp);
	}

	///////////////////////////////////////////////// Replacing values

	if (tIsAreaInsidePos == true)
	{
		if (basePos < posNewPoint)
		{
			replaceVectorRange(squarePoints, basePos, posNewPoint, pData, true);
		}
		else
		{
			replaceVectorRange(squarePoints, posNewPoint, basePos, pData, false);
		}
	}
	else
	{
		if (basePos < posNewPoint)
		{
			while (pData.size() > posNewPoint + 1)
			{
				pData.pop_back();
			}

			replaceVectorRange(squarePoints, 0, basePos, pData, false, false);
		}
		else
		{
			while (pData.size() > basePos + 1)
			{
				pData.pop_back();
			}

			replaceVectorRange(squarePoints, 0, posNewPoint, pData, true, false);
		}
	}
}

Line ImplantTools::GetSquareCornerFeatures(const Point& pBasePoint, const Point& pCenter, const Point& pSidePoint, const std::vector<Point>& pData, int& pDataSidePos, double pSlopeAngle)
{
	Point vector1 = pSidePoint - pCenter;
	Point vector2 = pBasePoint - pCenter;
	double angle = ImplantTools::getAngleBetweenVectorsDegree(vector1, vector2);

	Point newPoint = pBasePoint;

	Line myLine = Line::makeLineWithPoints(pCenter, pSidePoint);

	if (angle < pSlopeAngle - 0.5)
	{
		//applying the law of sines.
		double myAngle = 180. - pSlopeAngle;
		double dist1 = sqrt(vector1.dot(vector1));
		double dist2 = sqrt(vector2.dot(vector2));
		double angle2 = 180 - (angle + myAngle);
		double moveDist = (dist1 / sin(myAngle * PI / 180.)) * sin(angle2 * PI / 180.);

		vector2.normalice();
		if (moveDist < dist2)
		{
			newPoint = pCenter + moveDist * vector2;
		}
		else
		{
			newPoint = pBasePoint;
		}

		myLine.setDirectVector((newPoint - pSidePoint));
		myLine.setPoint(pCenter);
	}

	auto it1 = pData.begin();
	auto it2 = pData.end();
	double distance, secondDist = -1;

	int cont = 0;
	for (; it1 != it2; ++it1)
	{
		distance = myLine.getDistanceFromPoint(*it1);
		if (distance < secondDist || secondDist < 0)
		{
			secondDist = distance;
			newPoint = *it1;
			pDataSidePos = cont;
		}
		cont++;
	}

	myLine.setPoint(newPoint);
	return myLine;
}

void ImplantTools::squareCorner(const Point& pBasePoint, const Point& pCenter, const Point& pSecondPoint, std::vector<Point>& pData, double pSlopeAngle)
{
	Plane basePlane;
	Point vector1 = pBasePoint - pCenter;
	Point vector2 = pSecondPoint - pCenter;
	vector1.normalice();
	vector2.normalice();

	Point farBase = pBasePoint + 1000 * vector1;
	Point farSecond = pSecondPoint + 1000 * vector2;

	basePlane.init(vector1, farBase);
	Line myLine(vector1, farSecond);

	////////////////////////////////////////////Find extreme points

	auto it1 = pData.begin();
	auto it2 = pData.end();

	double distance, baseDist = -1, secondDist = -1;

	Point tBasePoint, tSecondPoint;

	int basePos, secondPos, cont = 0;

	for (; it1 != it2; it1++)
	{
		distance = basePlane.getDistanceFromPoint(*it1);

		if (distance < baseDist || baseDist == -1)
		{
			baseDist = distance;
			tBasePoint = *it1;
			basePos = cont;
		}

		distance = myLine.getDistanceFromPoint(*it1);

		if (distance < secondDist || secondDist == -1)
		{
			secondDist = distance;
			tSecondPoint = *it1;
			secondPos = cont;
		}

		cont++;
	}

	bool tIsAreaInsidePos;

	//////////////////////////////////////////
	/*
	looking for the area with the greatest number of points to avoid false positives with the areas of belonging between the planes
	*/

	int distInside, distBegin, distEnd;
	distInside = abs(basePos - secondPos);

	if (basePos < secondPos)
	{
		distBegin = basePos;
		distEnd = pData.size() - secondPos;
	}

	else
	{
		distBegin = secondPos;
		distEnd = pData.size() - basePos;
	}

	//////////////////////////////////// Detecting where the area to square is in the pData vector.

	Plane currentPlane, pBase, pSecond;
	currentPlane.init(pBasePoint, pCenter, pSecondPoint);
	pBase = currentPlane.getPerpendicularPlane(tSecondPoint, pCenter);
	pSecond = currentPlane.getPerpendicularPlane(tBasePoint, pCenter);

	pBase.reverseByPoint(farBase);
	pSecond.reverseByPoint(farSecond);

	if (distInside >= distBegin && distBegin >= distEnd)
	{
		int evalPos = (basePos + secondPos) / 2;

		if (pBase.eval(pData[evalPos]) > 0 && pSecond.eval(pData[evalPos]) > 0)
		{
			tIsAreaInsidePos = true;
		}
		else
		{
			tIsAreaInsidePos = false;
		}
	}
	else
	{
		int evalPos;
		int mid = (distEnd + distBegin) / 2;

		if (mid < distEnd)
		{
			evalPos = pData.size() - distEnd + mid;
		}
		else
		{
			evalPos = mid - distEnd;
		}

		if (pBase.eval(pData[evalPos]) > 0 && pSecond.eval(pData[evalPos]) > 0)
		{
			tIsAreaInsidePos = false;
		}
		else
		{
			tIsAreaInsidePos = true;
		}
	}

	/////////////////////////////////////////// Generate the new corner
	basePlane.movePlane(tBasePoint);
	myLine.setPoint(tSecondPoint);
	Point intercepTop = basePlane.getInterceptionLinePoint(myLine);
	Point lastTopPoint = tBasePoint;
	vector1 = tSecondPoint - intercepTop;
	vector2 = lastTopPoint - intercepTop;
	double angle = ImplantTools::getAngleBetweenVectorsDegree(vector1, vector2);

	Point newPoint = tSecondPoint;
	int posNewPoint = secondPos;

	if (angle < pSlopeAngle - 0.5)
	{
		//applying the law of sines.
		double myAngle = 180. - pSlopeAngle;
		double dist1 = sqrt(vector1.dot(vector1));
		double dist2 = sqrt(vector2.dot(vector2));
		double angle2 = 180 - (angle + myAngle);
		double moveDist = (dist1 / sin(myAngle * PI / 180.)) * sin(angle2 * PI / 180.);

		vector2.normalice();
		if (moveDist < dist2)
		{
			newPoint = intercepTop + moveDist * vector2;
		}
		else
		{
			newPoint = lastTopPoint;
		}

		myLine.setDirectVector((newPoint - tSecondPoint));
		myLine.setPoint(intercepTop);

		it1 = pData.begin();
		it2 = pData.end();
		secondDist = -1;

		cont = 0;
		for (; it1 != it2; ++it1)
		{
			distance = myLine.getDistanceFromPoint(*it1);
			if (distance < secondDist || secondDist < 0)
			{
				secondDist = distance;
				newPoint = *it1;
				posNewPoint = cont;
			}
			cont++;
		}

		myLine.setPoint(newPoint);
		intercepTop = basePlane.getInterceptionLinePoint(myLine);
	}

	std::vector<Point> squarePoints, tempAll;

	if (intercepTop != lastTopPoint)
	{
		for (double i = 0.2; i < 0.9; i += 0.2)
		{
			Point temp = lastTopPoint + i * (intercepTop - lastTopPoint);
			squarePoints.push_back(temp);
		}
		squarePoints.push_back(intercepTop);
	}

	for (double i = 0.2; i < 0.9; i += 0.2)
	{
		Point temp = intercepTop + i * (newPoint - intercepTop);
		squarePoints.push_back(temp);
	}

	///////////////////////////////////////////////// Replacing values

	if (tIsAreaInsidePos == true)
	{
		if (basePos < posNewPoint)
		{
			replaceVectorRange(squarePoints, basePos, posNewPoint, pData, true);
		}
		else
		{
			replaceVectorRange(squarePoints, posNewPoint, basePos, pData, false);
		}
	}
	else
	{
		if (basePos < posNewPoint)
		{
			while (pData.size() > posNewPoint + 1)
			{
				pData.pop_back();
			}

			replaceVectorRange(squarePoints, 0, basePos, pData, false);
		}
		else
		{
			while (pData.size() > basePos + 1)
			{
				pData.pop_back();
			}

			replaceVectorRange(squarePoints, 0, posNewPoint, pData, true);
		}
	}
}

Point ImplantTools::getHighestPointsOnTibia(const std::vector<Point>& pData, const Plane& pSplitLatPlane, Line& pTopLine, int& latPos, int& medPos)
{
	double distLatTemp, distMedTemp;
	int posLatTemp = -1, posMedTemp = -1;
	latPos = -1;
	medPos = -1;

	int tHullSize = pData.size();

	if (tHullSize == 0)
	{
		return Point();
	}

	Point hullCenter;
	Point highestPoint = pTopLine.getPoint();
	int iterations = 10;

	do
	{
		hullCenter = Point();
		double distLat = -1, distMed = -1;

		for (int i = 0; i < tHullSize; i++)
		{
			Point tPoint = pData[i];

			hullCenter = hullCenter + tPoint;

			if (pSplitLatPlane.eval(tPoint) > 0)
			{
				distLatTemp = pTopLine.getDistanceFromPoint(tPoint);

				if (distLat < 0 || distLat > distLatTemp)
				{
					distLat = distLatTemp;
					latPos = i;
				}
			}
			else if (pSplitLatPlane.eval(tPoint) < 0)
			{
				distMedTemp = pTopLine.getDistanceFromPoint(tPoint);

				if (distMed < 0 || distMed > distMedTemp)
				{
					distMed = distMedTemp;
					medPos = i;
				}
			}
		}

		if (latPos != posLatTemp && medPos != posMedTemp)
		{
			posLatTemp = latPos;
			posMedTemp = medPos;
			pTopLine = Line::makeLineWithPoints(pData[latPos], pData[medPos]);
			pTopLine.setPoint(highestPoint);
		}
		else
		{
			break;
		}
		iterations--;

	} while (iterations > 0);

	hullCenter = hullCenter / double(tHullSize);
	pTopLine = Line::makeLineWithPoints(pData[latPos], pData[medPos]);

	return hullCenter;
}

Plane ImplantTools::TransformPlane(const Plane& plane, const itk::Rigid3DTransform<>::Pointer pTransform)
{
	cv::Mat rotation = Rigid3DTransformToCVRotation(pTransform);
	cv::Mat translation = Rigid3DTransformToCVTranslation(pTransform);

	cv::Mat transformNormalVector = rotation * plane.getNormalVectorMat();
	cv::Mat transformPoint = (rotation * plane.getPointMat()) + translation;
	Plane transformPlane;
	transformPlane.init(Point(transformNormalVector), Point(transformPoint));
	return transformPlane;
}

std::vector<Point> ImplantTools::getSortPointVectorFill(const std::vector<Point>& pVector, double pDistance)
{
	if (pVector.size() < 2)
	{
		return pVector;
	}

	std::vector<Point> result;

	auto it1 = pVector.begin() + 1;
	auto it2 = pVector.end();

	result.push_back(*pVector.begin());
	int realSizePos = 0;

	for (; it1 != it2; ++it1)
	{
		Point a = result[realSizePos];
		Point b = *it1;

		double tDist = getDistanceBetweenPoints(a, b);

		if (tDist > pDistance)
		{
			Point vector = b - a;
			vector.normalice();
			int div = (tDist / pDistance) - 1;
			while (div > 0)
			{
				result.push_back((a + pDistance * vector));
				realSizePos++;
				div--;
			}
			result.push_back(b);
			realSizePos++;
		}
		else
		{
			result.push_back(b);
			realSizePos++;
		}

	}

	return result;
}

bool ImplantTools::areBothSetOfPointsSeparated(const std::vector<Point>& pPoints1, const std::vector<Point>& pPoints2, const cv::Mat& myRotationZ, float margin)
{
	if (pPoints1.size() == 0 || pPoints2.size() == 0)
	{
		return true;
	}

	ConvexHull hullObj1(pPoints1, myRotationZ);
	ConvexHull hullObj2(pPoints2, myRotationZ);

	std::vector<Point> convexHull1 = hullObj1.GetConvexHull();
	std::vector<Point> convexHull2 = hullObj2.GetConvexHull();

	bool result1 = hullObj1.areSomePointWithinConvexHull(convexHull2, margin);
	bool result2 = hullObj2.areSomePointWithinConvexHull(convexHull1, margin);
	bool areIntersection = result1 || result2;

	return !areIntersection;
}

std::vector<Point> ImplantTools::increaseVectorPoints(const std::vector<Point>& pPoints, int beginPos, int endPos, float distance)
{
	std::vector<Point> result;
	int tSize = pPoints.size();
	for (int i = 0; i < tSize; i++)
	{
		if ((i >= beginPos) && (i + 1 <= endPos) && (i + 1 < tSize))
		{
			result.push_back(pPoints[i]);
			float dist = ImplantTools::getDistanceBetweenPoints(pPoints[i], pPoints[i + 1]);
			if (dist > distance)
			{
				int diff = dist - distance;
				if (diff > 1)
				{
					float weight = 1. / float(diff);
					for (float j = 1; j < diff; j++)
					{
						Point vector = pPoints[i + 1] - pPoints[i];
						result.push_back(pPoints[i] + j * weight * vector);
					}
				}
			}
		}
		else
		{
			result.push_back(pPoints[i]);
		}
	}
	return result;
}

vtkSmartPointer<vtkPolyData> ImplantTools::getPolyLine(const std::vector<Point>& sortPoints)
{
	vtkNew<vtkPoints> points;
	for (int i = 0; i < sortPoints.size(); i++)
	{
		double pnt[3];
		pnt[0] = sortPoints[i].x;
		pnt[1] = sortPoints[i].y;
		pnt[2] = sortPoints[i].z;

		points->InsertNextPoint(pnt);
	}

	vtkNew<vtkCellArray> cells;

	for (unsigned int i = 0; i < points->GetNumberOfPoints() - 1; i++)
	{
		vtkNew<vtkLine> myLine;

		myLine->GetPointIds()->SetId(0, i);
		myLine->GetPointIds()->SetId(1, i + 1);

		cells->InsertNextCell(myLine);
	}

	vtkNew<vtkPolyData> polyData;
	polyData->SetPoints(points);
	polyData->SetLines(cells);
	return polyData;
}

void ImplantTools::show(const vtkSmartPointer<vtkPolyData> poly1, const vtkSmartPointer<vtkPolyData> poly2)
{
	vtkNew<vtkNamedColors> colors;

	vtkNew<vtkPolyDataMapper> contoursMapper;
	contoursMapper->SetInputData(poly1);
	contoursMapper->ScalarVisibilityOff();

	vtkNew<vtkActor> contoursActor;
	contoursActor->SetMapper(contoursMapper);
	contoursActor->GetProperty()->SetRepresentationToWireframe();
	contoursActor->GetProperty()->ShadingOff();
	contoursActor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

	vtkNew<vtkPolyDataMapper> surfaceMapper;
	surfaceMapper->SetInputData(poly2);
	surfaceMapper->ScalarVisibilityOff();

	vtkNew<vtkActor> surfaceActor;
	surfaceActor->SetMapper(surfaceMapper);
	surfaceActor->GetProperty()->SetRepresentationToWireframe();
	surfaceActor->GetProperty()->ShadingOff();
	surfaceActor->GetProperty()->SetColor(
		colors->GetColor3d("MistyRose").GetData());

	// Create two renderers side by side to show the contours and the surface
	// separately
	//
	std::cout << "Press 't' for trackball interaction" << std::endl;
	std::cout << "Press 'r' to reset the camera" << std::endl;
	std::cout << "Press 'w' for wireframe representation" << std::endl;
	std::cout << "Press 's' for surface representation" << std::endl;

	vtkNew<vtkRenderer> renderer1;
	renderer1->SetViewport(0., 0., 0.5, 1.);
	renderer1->SetBackground(colors->GetColor3d("CadetBlue").GetData());

	vtkNew<vtkRenderer> renderer2;
	renderer2->SetViewport(0.5, 0., 1., 1.);
	renderer2->SetBackground(colors->GetColor3d("BurlyWood").GetData());
	//renderer2->SetBackground(colors->GetColor3d("SlateGray").GetData());

	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->SetSize(800, 400);
	renderWindow->SetWindowName("ContoursToSurface");

	renderWindow->AddRenderer(renderer1);
	renderWindow->AddRenderer(renderer2);

	vtkNew<vtkRenderWindowInteractor> interactor;
	interactor->SetRenderWindow(renderWindow);

	renderer1->AddViewProp(surfaceActor);
	renderer2->AddViewProp(contoursActor);
	renderWindow->Render();

	interactor->Start();
}


void ImplantTools::show(vtkSmartPointer<vtkPolyData> poly, const std::vector<Point>& points, bool makePolyLine)
{
	vtkNew<vtkNamedColors> colors;

	vtkNew<vtkPolyDataMapper> contoursMapper;
	contoursMapper->SetInputData(poly);
	contoursMapper->ScalarVisibilityOff();

	vtkNew<vtkActor> polyActor;
	polyActor->SetMapper(contoursMapper);
	polyActor->GetProperty()->SetRepresentationToWireframe();
	polyActor->GetProperty()->ShadingOff();
	polyActor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

	std::vector<vtkSmartPointer<vtkActor>> pointsActor;

	if (makePolyLine == true)
	{
		vtkSmartPointer<vtkPolyData> polyLine = getPolyLine(points);

		vtkNew<vtkPolyDataMapper> polyMapper;
		polyMapper->SetInputData(polyLine);
		polyMapper->ScalarVisibilityOff();

		vtkNew<vtkActor> polyActor;
		polyActor->SetMapper(polyMapper);
		polyActor->GetProperty()->SetRepresentationToWireframe();
		polyActor->GetProperty()->ShadingOff();
		polyActor->GetProperty()->SetColor(colors->GetColor3d("blue").GetData());

		pointsActor.push_back(polyActor);
	}
	else
	{
		for (int i = 0; i < points.size(); i++)
		{
			double pnt[3];
			pnt[0] = points[i].x;
			pnt[1] = points[i].y;
			pnt[2] = points[i].z;

			vtkNew<vtkSphereSource> sphere;
			sphere->SetCenter(pnt);
			sphere->SetRadius(0.2);
			sphere->Update();

			vtkNew<vtkPolyDataMapper> sphereMapper;
			sphereMapper->SetInputData(sphere->GetOutput());
			sphereMapper->ScalarVisibilityOff();

			vtkNew<vtkActor> sphereActor;
			sphereActor->SetMapper(sphereMapper);
			sphereActor->GetProperty()->SetRepresentationToWireframe();
			sphereActor->GetProperty()->ShadingOff();
			sphereActor->GetProperty()->SetColor(colors->GetColor3d("blue").GetData());

			pointsActor.push_back(sphereActor);
		}
	}

	vtkNew<vtkRenderer> renderer;
	//renderer->SetViewport(0., 0., 0.5, 1.);
	renderer->SetBackground(colors->GetColor3d("CadetBlue").GetData());

	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->SetSize(800, 400);
	renderWindow->SetWindowName("Surface");

	renderWindow->AddRenderer(renderer);

	vtkNew<vtkRenderWindowInteractor> interactor;
	interactor->SetRenderWindow(renderWindow);

	renderer->AddActor(polyActor);

	for (int i = 0; i < pointsActor.size(); i++)
	{
		renderer->AddActor(pointsActor[i]);
	}

	renderWindow->Render();

	interactor->Start();
}

