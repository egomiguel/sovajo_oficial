#include "HipFemurRegistrationAnterior.hpp"
#include "RegistrationPrivate.hpp"
#include "TemplatePointsHip.hpp"
#include "LeastSquaresICP.hpp"
#include "RegistrationException.hpp"

using namespace THA::RIGISTRATION;

HipRegistrationFemurAnterior::HipRegistrationFemurAnterior(const vtkSmartPointer<vtkPolyData> pImage, const PointTypeITK& pAnteriorFemoralNeckCT, const PointTypeITK& pAnteriorDistalTrochanterCT, const PointTypeITK& pLateralTrochanterCT, RegisterSide pSide)
    :HipRegistrationFemur(pImage, pAnteriorFemoralNeckCT, pAnteriorDistalTrochanterCT, pLateralTrochanterCT, pSide)
{
}

bool HipRegistrationFemurAnterior::RegistrationLandmarksAnterior(const PointTypeITK& pAnteriorFemoralNeckCamera, const PointTypeITK& pAnteriorDistalTrochanterCamera, const PointTypeITK& pLateralTrochanterCamera, double& error)
{
	return HipRegistrationFemur::RegistrationLandmarks(pAnteriorFemoralNeckCamera, pAnteriorDistalTrochanterCamera, pLateralTrochanterCamera, error);
}

bool HipRegistrationFemurAnterior::MakeRegistrationAnterior(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pAnteriorFemoralNeckCamera, const PointTypeITK& pAnteriorDistalTrochanterCamera, const PointTypeITK& pLateralTrochanterCamera)
{
	return HipRegistrationFemur::MakeRegistration(pBonePoints, pAnteriorFemoralNeckCamera, pAnteriorDistalTrochanterCamera, pLateralTrochanterCamera);
}

std::vector<RegistrationPointsHip> HipRegistrationFemurAnterior::getRegistrationPointAnteriorLateral(std::vector<PointTypeITK>& pVerificationPoints, double& pError) const
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

std::vector<RegistrationPointsHip> HipRegistrationFemurAnterior::getRegistrationPointAnterior(std::vector<PointTypeITK>& pVerificationPoints, double& pError) const
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

cv::Mat HipRegistrationFemurAnterior::getTemplateAlignment(double& pError) const
{
	TemplateHipFemur templateObj;

	cv::Point3d targetCenter;
	cv::Point3d sourceCenter;

	std::vector<cv::Point3d> vectorSource, vectorTarget, myTemplatePoints;

	targetCenter = mCenter;
	vectorTarget = m_data->getAxisHipFemoral(mFemoralNeck, mDistalTrochanter, mLateralTrochanter);

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

