#include "HipPelvisRegistration.hpp"
#include "RegistrationPrivate.hpp"
#include "TemplatePointsHip.hpp"
#include "LeastSquaresICP.hpp"
#include "RegistrationException.hpp"

using namespace THA::RIGISTRATION;

PelvisRegistration::PelvisRegistration(const vtkSmartPointer<vtkPolyData> pImage, const PointTypeITK& pAnteriorAcetabulum, const PointTypeITK& pPosteriorAcetabulum, const PointTypeITK& pSuperiorAcetabulum, RegisterSide pSide)
	:Registration(pImage)
{

	vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
	implicitPolyDataDistance->SetInput(Registration::poly);

	/*double signedDistance1 = m_data->GetNearPoint(implicitPolyDataDistance, pAnteriorAcetabulum, anteriorAcetabulum);
	double signedDistance2 = m_data->GetNearPoint(implicitPolyDataDistance, pPosteriorAcetabulum, posteriorAcetabulum);
	double signedDistance3 = m_data->GetNearPoint(implicitPolyDataDistance, pSuperiorAcetabulum, superiorAcetabulum);

	if (signedDistance1 > 5 || signedDistance2 > 5 || signedDistance3 > 5)
	{
		throw RegistrationException("Check your CT landmarks. Distance from the bone model is too large.");
	}*/

	anteriorAcetabulum = pAnteriorAcetabulum;
	posteriorAcetabulum = pPosteriorAcetabulum;
	superiorAcetabulum = pSuperiorAcetabulum;

	cv::Point3d anteriorCV, posteriorCV, superiorCV;
	anteriorCV = m_data->itkPointToCV(anteriorAcetabulum);
	posteriorCV = m_data->itkPointToCV(posteriorAcetabulum);
	superiorCV = m_data->itkPointToCV(superiorAcetabulum);

	auto minCircle = m_data->getMinCircle(anteriorCV, posteriorCV, superiorCV);

	mCenter = minCircle.second;
	mRadius = minCircle.first;

	cv::Point3d refVectorAnterior = anteriorCV - mCenter;
	cv::Point3d refVectorSuperior = superiorCV - mCenter;
	cv::Point3d refVector = refVectorSuperior.cross(refVectorAnterior);
	mVectorNormal = refVector / sqrt(refVector.dot(refVector));

	mSide = pSide;

	/*
	cv::Point3d refA, refB;
	refA = mCenter + mRadius * mVectorNormal;
	refB = mCenter - mRadius * mVectorNormal;

	double pntA[3] = { refA.x, refA.y, refA.z };
	double pntB[3] = { refB.x, refB.y, refB.z };

	double myClosest[3];

	RegisterSide tempSise;
	if (abs(implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pntA, myClosest)) < abs(implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pntB, myClosest)))
	{
		tempSise = RegisterSide::RIGHT;
	}
	else
	{
		tempSise = RegisterSide::LEFT;
	}

	if (pSide == RegisterSide::UNKNOWN)
	{
		mSide = tempSise;
	}
	else
	{
		if (pSide != tempSise && pSide == RegisterSide::RIGHT)
		{
			throw RegistrationExceptionCode::REFERENCE_POINTS_DO_NOT_CORRESPOND_TO_PELVIS_RIGHT_SIDE;
		}
		else if (pSide != tempSise && pSide == RegisterSide::LEFT)
		{
			throw RegistrationExceptionCode::REFERENCE_POINTS_DO_NOT_CORRESPOND_TO_PELVIS_LEFT_SIDE;
		}
		else
		{
			mSide = pSide;
		}
	}*/
}

std::vector<RegistrationPointsHip> PelvisRegistration::getRegistrationPointPelvis(std::vector<PointTypeITK>& pVerificationPoints, double& pError) const
{
	vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
	implicitPolyDataDistance->SetInput(Registration::poly);

	std::vector<RegistrationPointsHip> result;

	TemplateHipPelvis templateObj;

	cv::Point3d targetSuperior, targetNormal, targetAxis, targetCenter;
	cv::Point3d sourceSuperior, sourceNormal, sourceAxis, sourceCenter;
	double tScale;
	std::vector<cv::Point3d> myTemplatePoints;

	targetSuperior = m_data->itkPointToCV(superiorAcetabulum) - mCenter;
	targetSuperior = targetSuperior / sqrt(targetSuperior.dot(targetSuperior));
	targetNormal = mVectorNormal;
	targetAxis = targetNormal.cross(targetSuperior);
	targetAxis = targetAxis / sqrt(targetAxis.dot(targetAxis));
	targetCenter = mCenter;

	std::vector<std::vector<cv::Point3d>> myOfficialPoints;
	std::vector<cv::Point3d> myVerificationPoints;

	if (mSide == RegisterSide::LEFT)
	{
		sourceSuperior = templateObj.vectorSuperiorLeft;
		sourceNormal = templateObj.vectorNormalLeft;
		sourceCenter = templateObj.centerLeft;
		sourceAxis = sourceNormal.cross(sourceSuperior);
		sourceAxis = sourceAxis / sqrt(sourceAxis.dot(sourceAxis));

		tScale = (mRadius / templateObj.sizeRadiusLeft);
		myTemplatePoints = templateObj.mTemplatePointsLeft;
		myOfficialPoints = templateObj.mTemplateOfficialPointsLeft;
		myVerificationPoints = templateObj.mTemplateVerificationPointsLeft;
	}
	else
	{
		sourceSuperior = templateObj.vectorSuperiorRight;
		sourceNormal = templateObj.vectorNormalRight;
		sourceCenter = templateObj.centerRight;
		sourceAxis = sourceNormal.cross(sourceSuperior);
		sourceAxis = sourceAxis / sqrt(sourceAxis.dot(sourceAxis));

		tScale = (mRadius / templateObj.sizeRadiusRight);
		myTemplatePoints = templateObj.mTemplatePointsRight;
		myOfficialPoints = templateObj.mTemplateOfficialPointsRight;
		myVerificationPoints = templateObj.mTemplateVerificationPointsRight;
	}

	std::vector<cv::Point3d> vectorTarget = { targetNormal, targetSuperior, targetAxis };
	std::vector<cv::Point3d> vectorSource = { sourceNormal, sourceSuperior, sourceAxis };

	cv::Mat data(7, 1, CV_64F);

	cv::Mat rotation = LeastSquaresICP::GetRotationAnglesXYZ(vectorSource, vectorTarget, data);
	cv::Mat translation = m_data->cvPointToMat(targetCenter) - (rotation * m_data->cvPointToMat(sourceCenter));

	data.at<double>(3, 0) = data.at<double>(0, 0);
	data.at<double>(4, 0) = data.at<double>(1, 0);
	data.at<double>(5, 0) = data.at<double>(2, 0);

	data.at<double>(0, 0) = translation.at<double>(0, 0);
	data.at<double>(1, 0) = translation.at<double>(1, 0);
	data.at<double>(2, 0) = translation.at<double>(2, 0);

	data.at<double>(6, 0) = tScale;

	LeastSquaresICP registerObj(myTemplatePoints);

	pError = registerObj.LeastSquaresScale(Registration::poly, data);

	////////////////////////////////////////////////////////////////////

	cv::Mat myTranslation(3, 1, CV_64F);
	myTranslation.at<double>(0, 0) = data.at<double>(0, 0);
	myTranslation.at<double>(1, 0) = data.at<double>(1, 0);
	myTranslation.at<double>(2, 0) = data.at<double>(2, 0);

	double angleX = data.at<double>(3, 0);
	double angleY = data.at<double>(4, 0);
	double angleZ = data.at<double>(5, 0);

	double scale = data.at<double>(6, 0);

	cv::Mat myRotation = registerObj.GetRotationMatrix(angleX, angleY, angleZ);

	for (int i = 0; i < myOfficialPoints.size(); i++)
	{
		std::vector<cv::Point3d> vectorTemp;
		for (int j = 0; j < myOfficialPoints[i].size(); j++)
		{
			cv::Mat transformPointMat = scale * (myRotation * (m_data->cvPointToMat(myOfficialPoints[i][j]))) + myTranslation;
			cv::Point3d transformPoint = cv::Point3d(transformPointMat);
			double myClosest[3];
			double pnt[3] = { transformPoint.x, transformPoint.y, transformPoint.z };
			implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);
			vectorTemp.push_back(cv::Point3d(myClosest[0], myClosest[1], myClosest[2]));
		}
		result.push_back(RegistrationPointsHip(vectorTemp));
	}

	for (int j = 0; j < myVerificationPoints.size(); j++)
	{
		cv::Mat transformPointMat = scale * (myRotation * (m_data->cvPointToMat(myVerificationPoints[j]))) + myTranslation;
		cv::Point3d transformPoint = cv::Point3d(transformPointMat);
		double myClosest[3];
		double pnt[3] = { transformPoint.x, transformPoint.y, transformPoint.z };
		implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);
		PointTypeITK tempPoint;
		tempPoint[0] = myClosest[0];
		tempPoint[1] = myClosest[1];
		tempPoint[2] = myClosest[2];
		pVerificationPoints.push_back(tempPoint);
	}

	return result;
}

bool PelvisRegistration::RegistrationLandmarks(const PointTypeITK& pAnteriorAcetabulumCamera, const PointTypeITK& pPosteriorAcetabulumCamera, const PointTypeITK& pSuperiorAcetabulumCamera, double& error)
{
	vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
	implicitPolyDataDistance->SetInput(Registration::poly);

	std::vector<PointTypeITK> source, target;

	cv::Point3d pAnteriorCamera = m_data->itkPointToCV(pAnteriorAcetabulumCamera);
	cv::Point3d pPosteriorCamera = m_data->itkPointToCV(pPosteriorAcetabulumCamera);
	cv::Point3d pSuperiorCamera = m_data->itkPointToCV(pSuperiorAcetabulumCamera);
	cv::Point3d pAnteriorCT = m_data->itkPointToCV(anteriorAcetabulum);
	cv::Point3d pPosteriorCT = m_data->itkPointToCV(posteriorAcetabulum);
	cv::Point3d pSuperiorCT = m_data->itkPointToCV(superiorAcetabulum);

	auto threeVectors = m_data->getAxisPelvisVectors(pAnteriorCamera, pPosteriorCamera, pSuperiorCamera, pAnteriorCT, pPosteriorCT, pSuperiorCT);
	auto minCircleSource = m_data->getMinCircle(pAnteriorCamera, pPosteriorCamera, pSuperiorCamera);

	Eigen::Matrix4d rigidTransform = m_data->getRigidTransform(threeVectors.first, threeVectors.second, minCircleSource.second, mCenter);

	PointTypeITK temp;

	double a = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pAnteriorAcetabulumCamera, temp);
	double b = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pPosteriorAcetabulumCamera, temp);
	double c = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pSuperiorAcetabulumCamera, temp);

	std::vector<double> distances = { a, b, c };
	std::sort(distances.rbegin(), distances.rend());

	error = distances[0];

	if (error > 7)
	{
		return false;
	}
	else
	{
		return true;
	}
}

bool PelvisRegistration::MakeRegistration(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pAnteriorAcetabulumCamera, const PointTypeITK& pPosteriorAcetabulumCamera, const PointTypeITK& pSuperiorAcetabulumCamera)
{
	std::vector<PointTypeITK> source, target;

	cv::Point3d pAnteriorCamera = m_data->itkPointToCV(pAnteriorAcetabulumCamera);
	cv::Point3d pPosteriorCamera = m_data->itkPointToCV(pPosteriorAcetabulumCamera);
	cv::Point3d pSuperiorCamera = m_data->itkPointToCV(pSuperiorAcetabulumCamera);
	cv::Point3d pAnteriorCT = m_data->itkPointToCV(anteriorAcetabulum);
	cv::Point3d pPosteriorCT = m_data->itkPointToCV(posteriorAcetabulum);
	cv::Point3d pSuperiorCT = m_data->itkPointToCV(superiorAcetabulum);

	auto threeVectors = m_data->getAxisPelvisVectors(pAnteriorCamera, pPosteriorCamera, pSuperiorCamera, pAnteriorCT, pPosteriorCT, pSuperiorCT);
	auto minCircleSource = m_data->getMinCircle(pAnteriorCamera, pPosteriorCamera, pSuperiorCamera);

	Eigen::Matrix4d rigidTransform = m_data->getRigidTransform(threeVectors.first, threeVectors.second, minCircleSource.second, mCenter);

	cv::Mat data(6, 1, CV_64F);

	Eigen::Matrix3d rotation(3, 3);

	rotation(0, 0) = rigidTransform(0, 0);
	rotation(1, 0) = rigidTransform(1, 0);
	rotation(2, 0) = rigidTransform(2, 0);

	rotation(0, 1) = rigidTransform(0, 1);
	rotation(1, 1) = rigidTransform(1, 1);
	rotation(2, 1) = rigidTransform(2, 1);

	rotation(0, 2) = rigidTransform(0, 2);
	rotation(1, 2) = rigidTransform(1, 2);
	rotation(2, 2) = rigidTransform(2, 2);

	Eigen::Vector3d rot = rotation.eulerAngles(0, 1, 2);

	data.at<double>(0, 0) = rigidTransform(0, 3);
	data.at<double>(1, 0) = rigidTransform(1, 3);
	data.at<double>(2, 0) = rigidTransform(2, 3);
	data.at<double>(3, 0) = rot(0);
	data.at<double>(4, 0) = rot(1);
	data.at<double>(5, 0) = rot(2);

	LeastSquaresICP myICP(pBonePoints);

	double error;

	error = myICP.LeastSquares(Registration::poly, data);

	Registration::MakeResult(data, error);

	return true;
}


/*
double PelvisRegistration::getAxisReference(cv::Point3d& pVectorOutSideNormal, cv::Point3d& pVectorSuperiorRing, cv::Point3d& pCenter)
{
	cv::Point3d anteriorCV, posteriorCV, superiorCV;
	anteriorCV = m_data->itkPointToCV(anteriorAcetabulum);
	posteriorCV = m_data->itkPointToCV(posteriorAcetabulum);
	superiorCV = m_data->itkPointToCV(superiorAcetabulum);

	std::vector<cv::Point3d> points = { anteriorCV, posteriorCV, superiorCV };

	RPlane helpPlane;
	helpPlane.init(anteriorCV, posteriorCV, superiorCV);

	double Z = 0;
	cv::Mat zRot = m_data->GetRotateZ(helpPlane.getNormalVector());
	std::vector<cv::Point2f> coplanar2d;

	auto it1 = points.begin();
	auto it2 = points.end();

	for (; it1 != it2; ++it1)
	{
		cv::Point3d pointTemp = m_data->TransformPoint(*it1, zRot);
		coplanar2d.push_back(cv::Point2f(pointTemp.x, pointTemp.y));
		Z += pointTemp.z;
	}

	Z = Z / double(points.size());

	cv::Point2f center;
	float radius;

	cv::minEnclosingCircle(coplanar2d, center, radius);

	cv::Point3d myCenter = m_data->TransformPoint(cv::Point3d(center.x, center.y, Z), zRot.inv());

	cv::Point3d tSuperiorVector = superiorCV - myCenter;
	pVectorSuperiorRing = tSuperiorVector / sqrt(tSuperiorVector.dot(tSuperiorVector));
	pCenter = myCenter;
	//////////////////////////////////

	cv::Point3d a = myCenter + radius * helpPlane.getNormalVector();
	cv::Point3d b = myCenter - radius * helpPlane.getNormalVector();

	double pntA[3] = { a.x, a.y, a.z };
	double pntB[3] = { b.x, b.y, b.z };

	vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
	implicitPolyDataDistance->SetInput(poly);
	double myClosest[3];
	cv::Point3d externalVector;

	if (abs(implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pntA, myClosest)) > abs(implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pntB, myClosest)))
	{
		helpPlane.reverseNormal(a, true);
	}
	else
	{
		helpPlane.reverseNormal(b, true);
	}
	pVectorOutSideNormal = helpPlane.getNormalVector();
	return radius;
}

std::vector<cv::Point3d> PelvisRegistration::getAcetabularPoints()
{
	int externalPoints = 15;

	cv::Point3d anteriorCV, posteriorCV, superiorCV;
	anteriorCV = m_data->itkPointToCV(anteriorAcetabulum);
	posteriorCV = m_data->itkPointToCV(posteriorAcetabulum);
	superiorCV = m_data->itkPointToCV(superiorAcetabulum);

	std::vector<cv::Point3d> points = { anteriorCV, posteriorCV, superiorCV };

	RPlane helpPlane;
	helpPlane.init(anteriorCV, posteriorCV, superiorCV);

	double Z = 0;
	cv::Mat zRot = m_data->GetRotateZ(helpPlane.getNormalVector());
	std::vector<cv::Point2f> coplanar2d;

	auto it1 = points.begin();
	auto it2 = points.end();

	for (; it1 != it2; ++it1)
	{
		cv::Point3d pointTemp = m_data->TransformPoint(*it1, zRot);
		coplanar2d.push_back(cv::Point2f(pointTemp.x, pointTemp.y));
		Z += pointTemp.z;
	}

	Z = Z / double(points.size());

	cv::Point2f center;
	float radius;

	cv::minEnclosingCircle(coplanar2d, center, radius);

	cv::Point3d myCenter = m_data->TransformPoint(cv::Point3d(center.x, center.y, Z), zRot.inv());

	cv::Point3d a = (anteriorCV + superiorCV) / 2.0;
	cv::Point3d b = (posteriorCV + superiorCV) / 2.0;

	RLine myLine = RLine::makeLineWithPoints(a, b);
	std::pair<cv::Point3d, cv::Point3d> interceptionPair;

	myLine.getInterceptionSphere(myCenter, radius, interceptionPair);

	a = interceptionPair.first;
	b = interceptionPair.second;

	cv::Point3d vectorA = a - myCenter;
	cv::Point3d vectorB = b - myCenter;
	cv::Point3d superiorVector = superiorCV - myCenter;

	RPlane splitPlane = helpPlane.getPerpendicularPlane(a, b);
	if (splitPlane.eval(superiorCV) < 0)
	{
		splitPlane.reverse();
	}

	double arc = RLine::getAngleBetweenVectors(vectorA, superiorVector) + RLine::getAngleBetweenVectors(vectorB, superiorVector);
	double step = arc / double(externalPoints - 1);

	cv::Mat rotMat = RLine::getRotateMatrix(helpPlane.getNormalVector(), step);
	cv::Point3d rotVector = m_data->TransformPoint(vectorA, rotMat);
	rotVector = rotVector / sqrt(rotVector.dot(rotVector));
	cv::Point3d pointTemp = myCenter + radius * rotVector;

	if (splitPlane.eval(pointTemp) < 0)
	{
		rotMat = RLine::getRotateMatrix(helpPlane.getNormalVector(), -step);
	}

	cv::Point3d moveVector = vectorA;
	std::vector<cv::Point3d> results;
	results.push_back(a);

	for (int i = 0; i < externalPoints - 2; i++)
	{
		rotVector = m_data->TransformPoint(moveVector, rotMat);
		rotVector = rotVector / sqrt(rotVector.dot(rotVector));
		pointTemp = myCenter + radius * rotVector;
		moveVector = rotVector;
		results.push_back(pointTemp);
	}
	results.push_back(b);

	a = myCenter + radius * helpPlane.getNormalVector();
	b = myCenter - radius * helpPlane.getNormalVector();

	double pntA[3] = {a.x, a.y, a.z};
	double pntB[3] = {b.x, b.y, b.z };

	vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
	implicitPolyDataDistance->SetInput(poly);
	double myClosest[3];
	cv::Point3d externalVector;

	if (implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pntA, myClosest) > implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pntB, myClosest))
	{
		helpPlane.reverseNormal(a, true);
	}
	else
	{
		helpPlane.reverseNormal(b, true);
	}

	it1 = results.begin();
	it2 = results.end();

	std::vector<cv::Point3d> borderPoints;

	cv::Point3d dowPoint = myCenter + 0.9 * radius * helpPlane.getNormalVector();

	int myCont = results.size() - 7;

	for (; it1 != it2; ++it1)
	{
		//borderPoints.push_back(*it1);

		pointTemp = (*it1) + radius * helpPlane.getNormalVector();

		RLine myLine = RLine::makeLineWithPoints(dowPoint, pointTemp);
		cv::Point3d vectorTemp = myCenter - (*it1);
		vectorTemp = vectorTemp / sqrt(vectorTemp.dot(vectorTemp));
		a = (*it1) + 5.0 * vectorTemp;
		b = (*it1) - 5.0 * vectorTemp;
		RPlane planeA, planeB, planeTemp;
		planeA.init(vectorTemp, a);
		planeB.init(vectorTemp, b);
		if (planeA.eval(*it1) < 0)
		{
			planeA.reverse();
		}

		if (planeB.eval(*it1) < 0)
		{
			planeB.reverse();
		}

		planeTemp = helpPlane.getPerpendicularPlane(dowPoint, pointTemp);

		vtkSmartPointer<vtkPolyData> contour = getContour(poly, planeTemp.getNormalVector(), planeTemp.getPoint());
		vtkSmartPointer<vtkPoints> pointsVTK = contour->GetPoints();
		int tSize = pointsVTK->GetNumberOfPoints();
		double distance = -1;
		cv::Point3d nearPoint;
		bool findPoint = false;
		vtkIdType nearPointId;
		for (int i = 0; i < tSize; i++)
		{
			double pnt[3];
			pointsVTK->GetPoint(i, pnt);
			if (planeA.eval(pnt) > 0 && planeB.eval(pnt) > 0)
			{
				findPoint = true;
				double distanceTemp = myLine.getDistanceFromPoint(pnt);
				if (distance < 0 || distance > distanceTemp)
				{
					distance = distanceTemp;
					nearPoint = cv::Point3d(pnt[0], pnt[1], pnt[2]);
					nearPointId = i;
				}
			}
		}

		if (findPoint == true)
		{
			std::vector<cv::Point3d> restPoints = m_data->getPointsAtDistanceVTK(nearPointId, 6, contour);

			double distA = m_data->GetDistance(restPoints[0], myCenter);
			double distB = m_data->GetDistance(restPoints[1], myCenter);

			if (distA > distB)
			{
				borderPoints.push_back(restPoints[0]);
			}
			else
			{
				borderPoints.push_back(restPoints[1]);
			}

			//for (int i = 0; i < restPoints.size(); i++)
			//{
			   // borderPoints.push_back(restPoints[i]);
			//}
			//borderPoints.push_back(nearPoint);
		}

		myCont--;
		if (myCont == 0)
		{
			break;
		}

	}

	return borderPoints;
}
*/