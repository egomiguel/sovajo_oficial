#include "HipFemurRegistration.hpp"
#include "RegistrationPrivate.hpp"
#include "TemplatePointsHip.hpp"
#include "LeastSquaresICP.hpp"
#include "RegistrationException.hpp"

using namespace THA::RIGISTRATION;

HipRegistrationFemur::HipRegistrationFemur(const vtkSmartPointer<vtkPolyData> pImage, const PointTypeITK& pAnteriorFemoralNeckCT, const PointTypeITK& pAnteriorDistalTrochanterCT, const PointTypeITK& pLateralTrochanterCT, RegisterSide pSide)
    :Registration(pImage)
{
    /*vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    implicitPolyDataDistance->SetInput(Registration::poly);

    double signedDistance1 = m_data->GetNearPoint(implicitPolyDataDistance, pAnteriorFemoralNeckCT, mAnteriorFemoralNeck);
    double signedDistance2 = m_data->GetNearPoint(implicitPolyDataDistance, pAnteriorDistalTrochanterCT, mAnteriorDistalTrochanter);
    double signedDistance3 = m_data->GetNearPoint(implicitPolyDataDistance, pLateralTrochanterCT, mLateralTrochanter);

    if (signedDistance1 > 5 || signedDistance2 > 5 || signedDistance3 > 5)
    {
        throw RegistrationException("Check your CT landmarks. Distance from the bone model is too large.");
    }*/

	mAnteriorFemoralNeck = pAnteriorFemoralNeckCT;
	mAnteriorDistalTrochanter = pAnteriorDistalTrochanterCT;
	mLateralTrochanter = pLateralTrochanterCT;

    mSide = pSide;

    auto minCircle = m_data->getMinCircle(mAnteriorFemoralNeck, mAnteriorDistalTrochanter, mLateralTrochanter);

    mCenter = minCircle.second;
}


bool HipRegistrationFemur::RegistrationLandmarks(const PointTypeITK& pAnteriorFemoralNeckCamera, const PointTypeITK& pAnteriorDistalTrochanterCamera, const PointTypeITK& pLateralTrochanterCamera, double& error)
{
    vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    implicitPolyDataDistance->SetInput(Registration::poly);

    auto source = m_data->getAxisHipFemoral(pAnteriorFemoralNeckCamera, pAnteriorDistalTrochanterCamera, pLateralTrochanterCamera);

    auto target = m_data->getAxisHipFemoral(mAnteriorFemoralNeck, mAnteriorDistalTrochanter, mLateralTrochanter);
    
    auto minCircleSource = m_data->getMinCircle(pAnteriorFemoralNeckCamera, pAnteriorDistalTrochanterCamera, pLateralTrochanterCamera);

    Eigen::Matrix4d rigidTransform = m_data->getRigidTransform(source, target, minCircleSource.second, mCenter);
    
    PointTypeITK temp;

    double a = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pAnteriorFemoralNeckCamera, temp);
    double b = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pAnteriorDistalTrochanterCamera, temp);
    double c = m_data->TransformPointDistance(implicitPolyDataDistance, rigidTransform, pLateralTrochanterCamera, temp);

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

bool HipRegistrationFemur::MakeRegistration(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pAnteriorFemoralNeckCamera, const PointTypeITK& pAnteriorDistalTrochanterCamera, const PointTypeITK& pLateralTrochanterCamera)
{
    auto source = m_data->getAxisHipFemoral(pAnteriorFemoralNeckCamera, pAnteriorDistalTrochanterCamera, pLateralTrochanterCamera);

    auto target = m_data->getAxisHipFemoral(mAnteriorFemoralNeck, mAnteriorDistalTrochanter, mLateralTrochanter);

    auto minCircleSource = m_data->getMinCircle(pAnteriorFemoralNeckCamera, pAnteriorDistalTrochanterCamera, pLateralTrochanterCamera);

    Eigen::Matrix4d rigidTransform = m_data->getRigidTransform(source, target, minCircleSource.second, mCenter);
    
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

std::vector<RegistrationPointsHip> HipRegistrationFemur::getRegistrationPointPosterolateral(std::vector<PointTypeITK>& pVerificationPoints, double& pError) const
{
    std::vector<RegistrationPointsHip> result;

    vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
    implicitPolyDataDistance->SetInput(Registration::poly);

    TemplateHipFemurOfficial templateOffcicial;

    std::vector<std::vector<cv::Point3d>> myOfficialPoints;
	std::vector<cv::Point3d> myVerification;

    if (mSide == RegisterSide::LEFT)
    {
        myOfficialPoints = templateOffcicial.mPosteroLateralLeft;
		myVerification = templateOffcicial.mVerificationPosteroLateralLeft;
    }
    else
    {
        myOfficialPoints = templateOffcicial.mPosteroLateralRight;
		myVerification = templateOffcicial.mVerificationPosteroLateralRight;
    }

    cv::Mat data = getTemplateAlignment(pError);

    cv::Mat myTranslation(3, 1, CV_64F);
    myTranslation.at<double>(0, 0) = data.at<double>(0, 0);
    myTranslation.at<double>(1, 0) = data.at<double>(1, 0);
    myTranslation.at<double>(2, 0) = data.at<double>(2, 0);

    double angleX = data.at<double>(3, 0);
    double angleY = data.at<double>(4, 0);
    double angleZ = data.at<double>(5, 0);

    double scale = data.at<double>(6, 0);

    std::vector<PointTypeITK> temp;
    LeastSquaresICP registerObj(temp);
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

	for (int j = 0; j < myVerification.size(); j++)
	{
		cv::Mat transformPointMat = scale * (myRotation * (m_data->cvPointToMat(myVerification[j]))) + myTranslation;
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

std::vector<RegistrationPointsHip> HipRegistrationFemur::getRegistrationPointAnterolateral(std::vector<PointTypeITK>& pVerificationPoints, double& pError) const
{
    std::vector<RegistrationPointsHip> result;

    vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
    implicitPolyDataDistance->SetInput(Registration::poly);

    TemplateHipFemurOfficial templateOffcicial;

    std::vector<std::vector<cv::Point3d>> myOfficialPoints;
	std::vector<cv::Point3d> myVerification;

    if (mSide == RegisterSide::LEFT)
    {
        myOfficialPoints = templateOffcicial.mAnteroLateralLeft;
		myVerification = templateOffcicial.mVerificationAnteroLateralLeft;
    }
    else
    {
        myOfficialPoints = templateOffcicial.mAnteroLateralRight;
		myVerification = templateOffcicial.mVerificationAnteroLateralRight;
    }

    cv::Mat data = getTemplateAlignment(pError);

    cv::Mat myTranslation(3, 1, CV_64F);
    myTranslation.at<double>(0, 0) = data.at<double>(0, 0);
    myTranslation.at<double>(1, 0) = data.at<double>(1, 0);
    myTranslation.at<double>(2, 0) = data.at<double>(2, 0);

    double angleX = data.at<double>(3, 0);
    double angleY = data.at<double>(4, 0);
    double angleZ = data.at<double>(5, 0);

    double scale = data.at<double>(6, 0);

    std::vector<PointTypeITK> temp;
    LeastSquaresICP registerObj(temp);
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

	for (int j = 0; j < myVerification.size(); j++)
	{
		cv::Mat transformPointMat = scale * (myRotation * (m_data->cvPointToMat(myVerification[j]))) + myTranslation;
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

std::vector<RegistrationPointsHip> HipRegistrationFemur::getRegistrationPointAnterior(std::vector<PointTypeITK>& pVerificationPoints, double& pError) const
{
    std::vector<RegistrationPointsHip> result;

    vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
    implicitPolyDataDistance->SetInput(Registration::poly);

    TemplateHipFemurOfficial templateOffcicial;

    std::vector<std::vector<cv::Point3d>> myOfficialPoints;
	std::vector<cv::Point3d> myVerification;

    if (mSide == RegisterSide::LEFT)
    {
        myOfficialPoints = templateOffcicial.mAnteriorLeft;
		myVerification = templateOffcicial.mVerificationAnteriorLeft;
    }
    else
    {
        myOfficialPoints = templateOffcicial.mAnteriorRight;
		myVerification = templateOffcicial.mVerificationAnteriorRight;
    }

    cv::Mat data = getTemplateAlignment(pError);

    cv::Mat myTranslation(3, 1, CV_64F);
    myTranslation.at<double>(0, 0) = data.at<double>(0, 0);
    myTranslation.at<double>(1, 0) = data.at<double>(1, 0);
    myTranslation.at<double>(2, 0) = data.at<double>(2, 0);

    double angleX = data.at<double>(3, 0);
    double angleY = data.at<double>(4, 0);
    double angleZ = data.at<double>(5, 0);

    double scale = data.at<double>(6, 0);

    std::vector<PointTypeITK> temp;
    LeastSquaresICP registerObj(temp);
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

	for (int j = 0; j < myVerification.size(); j++)
	{
		cv::Mat transformPointMat = scale * (myRotation * (m_data->cvPointToMat(myVerification[j]))) + myTranslation;
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

cv::Mat HipRegistrationFemur::getTemplateAlignment(double& pError) const
{
    TemplateHipFemur templateObj;

    cv::Point3d targetCenter;
    cv::Point3d sourceCenter;

    std::vector<cv::Point3d> vectorSource, vectorTarget, myTemplatePoints;

    targetCenter = mCenter;
    vectorTarget = m_data->getAxisHipFemoral(mAnteriorFemoralNeck, mAnteriorDistalTrochanter, mLateralTrochanter);

    if (mSide == RegisterSide::LEFT)
    {
        auto minCircleSource = m_data->getMinCircle(templateObj.mAnteriorNeckLeft, templateObj.mAnteriorDistalLeft, templateObj.mLateralTrochanterLeft);
        sourceCenter = minCircleSource.second;
        vectorSource = m_data->getAxisHipFemoral(templateObj.mAnteriorNeckLeft, templateObj.mAnteriorDistalLeft, templateObj.mLateralTrochanterLeft);

        myTemplatePoints = templateObj.mTemplatePointsLeft;
    }
    else
    {
        auto minCircleSource = m_data->getMinCircle(templateObj.mAnteriorNeckRight, templateObj.mAnteriorDistalRight, templateObj.mLateralTrochanterRight);
        sourceCenter = minCircleSource.second;
        vectorSource = m_data->getAxisHipFemoral(templateObj.mAnteriorNeckRight, templateObj.mAnteriorDistalRight, templateObj.mLateralTrochanterRight);

        myTemplatePoints = templateObj.mTemplatePointsRight;
    }

    cv::Mat data(7, 1, CV_64F);

    cv::Mat rotation = LeastSquaresICP::GetRotationAnglesXYZ(vectorSource, vectorTarget, data);
    cv::Mat translation = m_data->cvPointToMat(targetCenter) - (rotation * m_data->cvPointToMat(sourceCenter));

    data.at<double>(3, 0) = data.at<double>(0, 0);
    data.at<double>(4, 0) = data.at<double>(1, 0);
    data.at<double>(5, 0) = data.at<double>(2, 0);

    data.at<double>(0, 0) = translation.at<double>(0, 0);
    data.at<double>(1, 0) = translation.at<double>(1, 0);
    data.at<double>(2, 0) = translation.at<double>(2, 0);

    data.at<double>(6, 0) = 1.0;

    LeastSquaresICP registerObj(myTemplatePoints);

    pError = registerObj.LeastSquaresScale(Registration::poly, data);

    return data;
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