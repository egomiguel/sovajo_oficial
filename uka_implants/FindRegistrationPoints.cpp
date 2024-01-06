//#include "RegistrationPointsVTK.hpp"
#include "FindRegistrationPoints.hpp"
#include "LeastSquaresScaleICP.hpp"
#include "TemplatePoints.hpp"
#include "ImplantTools.hpp"

using namespace UKA::IMPLANTS;

inline std::vector<PointTypeITK> PointToITKVector(const std::vector<Point>& pData)
{
    std::vector<PointTypeITK> result;
    auto it1 = pData.begin();
    auto it2 = pData.end();

    for (; it1 != it2; ++it1)
    {
        PointTypeITK temp;
        temp[0] = (*it1).x;
        temp[1] = (*it1).y;
        temp[2] = (*it1).z;
        result.push_back(temp);
    }
    return result;
}

RegistrationPoints::RegistrationPoints(std::vector<PointTypeITK> pPoints)
{
    points = pPoints;
}

RegistrationPoints::RegistrationPoints(std::vector<Point> pPoints)
{
    points = PointToITKVector(pPoints);
}

FindRegistrationPoints::FindRegistrationPoints(const Knee& pKnee)
{
    mKnee = pKnee;
}

FindRegistrationPoints::~FindRegistrationPoints()
{
}

std::vector<RegistrationPoints> FindRegistrationPoints::GetRegistrationPointsFemur(std::vector<Point>& pCheckPoints, double& pError)
{
    vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
    implicitPolyDataDistance->SetInput(mKnee.GetFemurPoly());

    std::vector<RegistrationPoints> result;

    TemplateFemur templateObj;

    Point targetAP, targetTEA, targetAxis, targetCenter;
    Point sourceAP, sourceTEA, sourceAxis, sourceCenter;
    double tScale;

    targetAP = mKnee.getFemurDirectVectorAP();
    targetTEA = mKnee.getFemurVectorLateralTEA();
    targetAxis = targetAP.cross(targetTEA);
    targetAxis.normalice();
    targetCenter = (mKnee.getLateralEpicondyle() + mKnee.getMedialEpicondyle()) / 2.0;

    Point diff = mKnee.getLateralEpicondyle() - mKnee.getMedialEpicondyle();
    double targetSize = sqrt(diff.dot(diff));

    std::vector<std::vector<Point>> myOfficialPoints;
    std::vector<cv::Point3d> myTemplatePoints;
    std::vector<Point> myCheckPoints;

    if (mKnee.getIsRight() == false)
    {
        sourceAP = templateObj.vectorLeftAP;
        sourceTEA = templateObj.vectorLeftTEA;
        sourceCenter = templateObj.centerLeft;
        tScale = (targetSize / templateObj.sizeLeft);
        myTemplatePoints = templateObj.mTemplateLeft;

		if (mKnee.getSurgerySide() == SurgerySideEnum::KLateral)
		{
			myOfficialPoints = templateObj.mOfficialLeftLateralPoints;
			myCheckPoints = templateObj.mLeftLateralCheckPoints;
		}
		else
		{
			myOfficialPoints = templateObj.mOfficialLeftMedialPoints;
			myCheckPoints = templateObj.mLeftMedialCheckPoints;
		}
    }
    else
    {
        sourceAP = templateObj.vectorRightAP;
        sourceTEA = templateObj.vectorRightTEA;
        sourceCenter = templateObj.centerRight;
        tScale = (targetSize / templateObj.sizeLeft);
        myTemplatePoints = templateObj.mTemplateRight;

		if (mKnee.getSurgerySide() == SurgerySideEnum::KLateral)
		{
			myOfficialPoints = templateObj.mOfficialRightLateralPoints;
			myCheckPoints = templateObj.mRightLateralCheckPoints;
		}
		else
		{
			myOfficialPoints = templateObj.mOfficialRightMedialPoints;
			myCheckPoints = templateObj.mRightMedialCheckPoints;
		}
    }

    if (tScale < 1.0)
    {
        tScale = 1.0;
    }

    sourceAxis = sourceAP.cross(sourceTEA);
    sourceAxis.normalice();

    std::vector<cv::Point3d> vectorTarget = { targetAP, targetTEA, targetAxis };
    std::vector<cv::Point3d> vectorSource = { sourceAP, sourceTEA, sourceAxis };

    cv::Mat data(7, 1, CV_64F);

    cv::Mat rotation = LeastSquaresScaleICP::GetRotationAnglesXYZ(vectorSource, vectorTarget, data);
    cv::Mat translation = targetCenter.ToMatPoint() - (rotation * sourceCenter.ToMatPoint());

    data.at<double>(3, 0) = data.at<double>(0, 0);
    data.at<double>(4, 0) = data.at<double>(1, 0);
    data.at<double>(5, 0) = data.at<double>(2, 0);

    data.at<double>(0, 0) = translation.at<double>(0, 0);
    data.at<double>(1, 0) = translation.at<double>(1, 0);
    data.at<double>(2, 0) = translation.at<double>(2, 0);

    data.at<double>(6, 0) = tScale;

    LeastSquaresScaleICP registerObj(myTemplatePoints);

    pError = registerObj.LeastSquaresScale(mKnee.GetFemurPoly(), data);

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
        std::vector<Point> vectorTemp;
        for (int j = 0; j < myOfficialPoints[i].size(); j++)
        {
            cv::Mat transformPointMat = scale * (myRotation * (myOfficialPoints[i][j]).ToMatPoint()) + myTranslation;
            Point transformPoint = Point(transformPointMat);
            double myClosest[3];
            double pnt[3] = { transformPoint.x, transformPoint.y, transformPoint.z };
            implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);
            vectorTemp.push_back(Point(myClosest[0], myClosest[1], myClosest[2]));
        }
        result.push_back(RegistrationPoints(vectorTemp));
    }

    for (int i = 0; i < myCheckPoints.size(); i++)
    {
        cv::Mat transformPointMat = scale * (myRotation * (myCheckPoints[i]).ToMatPoint()) + myTranslation;
        Point transformPoint = Point(transformPointMat);
        double myClosest[3];
        double pnt[3] = { transformPoint.x, transformPoint.y, transformPoint.z };
        implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);
        pCheckPoints.push_back(Point(myClosest[0], myClosest[1], myClosest[2]));
    }

    return result;
}

std::vector<RegistrationPoints> FindRegistrationPoints::GetRegistrationPointsTibia(std::vector<Point>& pCheckPoints, double& pError)
{
    vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
    implicitPolyDataDistance->SetInput(mKnee.GetTibiaPoly());

    std::vector<RegistrationPoints> result;

    TemplateTibia templateObj;

    Point targetAP, targetTEA, targetAxis, targetCenter;
    Point sourceAP, sourceTEA, sourceAxis, sourceCenter;
    double tScale;
    std::vector<cv::Point3d> myTemplatePoints;

    std::pair<double, double> axisSize = mKnee.getTibiaAutomaticAxis(targetAP, targetTEA, targetCenter);
    targetAxis = targetAP.cross(targetTEA);
    targetAxis.normalice();

    std::vector<std::vector<Point>> myOfficialPoints;
    //std::vector<Point> myInterior, myExterior;
    std::vector<Point> myCheckPoints;

    if (mKnee.getIsRight() == false)
    {
        sourceAP = templateObj.vectorLeftAP;
        sourceTEA = templateObj.vectorLeftTEA;
        sourceCenter = templateObj.centerLeft;
        tScale = ((axisSize.first / templateObj.sizeLeftAP) + (axisSize.second / templateObj.sizeLeftTEA)) / 2.0;
        myTemplatePoints = templateObj.mTemplateLeft;

		if (mKnee.getSurgerySide() == SurgerySideEnum::KLateral)
		{
			myOfficialPoints = templateObj.mOfficialLeftLateralPoints;
			myCheckPoints = templateObj.mLeftLateralCheckPoints;
		}
		else
		{
			myOfficialPoints = templateObj.mOfficialLeftMedialPoints;
			myCheckPoints = templateObj.mLeftMedialCheckPoints;
		}

        //myInterior = templateObj.mOfficialLeftPointsInterior;
        //myExterior = templateObj.mOfficialLeftPointsExterior;
        
    }
    else
    {
        sourceAP = templateObj.vectorRightAP;
        sourceTEA = templateObj.vectorRightTEA;
        sourceCenter = templateObj.centerRight;
        tScale = ((axisSize.first / templateObj.sizeRightAp) + (axisSize.second / templateObj.sizeRightTea)) / 2.0;
        myTemplatePoints = templateObj.mTemplateRight;

		if (mKnee.getSurgerySide() == SurgerySideEnum::KLateral)
		{
			myOfficialPoints = templateObj.mOfficialRightLateralPoints;
			myCheckPoints = templateObj.mRightLateralCheckPoints;
		}
		else
		{
			myOfficialPoints = templateObj.mOfficialRightMedialPoints;
			myCheckPoints = templateObj.mRightMedialCheckPoints;
		}

        //myInterior = templateObj.mOfficialRightPointsInterior;
        //myExterior = templateObj.mOfficialRightPointsExterior;
        
    }

    if (tScale < 1.0)
    {
        tScale = 1.0;
    }

    sourceAxis = sourceAP.cross(sourceTEA);
    sourceAxis.normalice();

    std::vector<cv::Point3d> vectorTarget = { targetAP, targetTEA, targetAxis };
    std::vector<cv::Point3d> vectorSource = { sourceAP, sourceTEA, sourceAxis };

    cv::Mat data(7, 1, CV_64F);

    cv::Mat rotation = LeastSquaresScaleICP::GetRotationAnglesXYZ(vectorSource, vectorTarget, data);
    cv::Mat translation = targetCenter.ToMatPoint() - (rotation * sourceCenter.ToMatPoint());

    data.at<double>(3, 0) = data.at<double>(0, 0);
    data.at<double>(4, 0) = data.at<double>(1, 0);
    data.at<double>(5, 0) = data.at<double>(2, 0);

    data.at<double>(0, 0) = translation.at<double>(0, 0);
    data.at<double>(1, 0) = translation.at<double>(1, 0);
    data.at<double>(2, 0) = translation.at<double>(2, 0);

    data.at<double>(6, 0) = tScale;

    LeastSquaresScaleICP registerObj(myTemplatePoints);

    pError = registerObj.LeastSquaresScale(mKnee.GetTibiaPoly(), data);

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
        std::vector<Point> vectorTemp;
        for (int j = 0; j < myOfficialPoints[i].size(); j++)
        {
            cv::Mat transformPointMat = scale * (myRotation * (myOfficialPoints[i][j]).ToMatPoint()) + myTranslation;
            Point transformPoint = Point(transformPointMat);
            double myClosest[3];
            double pnt[3] = { transformPoint.x, transformPoint.y, transformPoint.z };
            implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);
            vectorTemp.push_back(Point(myClosest[0], myClosest[1], myClosest[2]));
        }
        result.push_back(RegistrationPoints(vectorTemp));
    }

	/*
    std::vector<Point> finalPoints;

    for (int i = 0; i < myInterior.size(); i++)
    {
        cv::Mat transformInterior = scale * (myRotation * (myInterior[i]).ToMatPoint()) + myTranslation;
        cv::Mat transformExterior = scale * (myRotation * (myExterior[i]).ToMatPoint()) + myTranslation;

        Point interiorPoint = Point(transformInterior);
        Point exteriorPoint = Point(transformExterior);
        Point onSurface;

        double tError = ImplantTools::GetInterceptionWithLine(implicitPolyDataDistance, exteriorPoint, interiorPoint, onSurface);

        if (tError > 1)
        {
            finalPoints.push_back(transformInterior);
        }
        else
        {
            finalPoints.push_back(onSurface);
        }
    }

    result.push_back(RegistrationPoints(finalPoints));
	*/

    for (int i = 0; i < myCheckPoints.size(); i++)
    {
        cv::Mat transformPointMat = scale * (myRotation * (myCheckPoints[i]).ToMatPoint()) + myTranslation;
        Point transformPoint = Point(transformPointMat);
        double myClosest[3];
        double pnt[3] = { transformPoint.x, transformPoint.y, transformPoint.z };
        implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);
        pCheckPoints.push_back(Point(myClosest[0], myClosest[1], myClosest[2]));
    }

    return result;
}
