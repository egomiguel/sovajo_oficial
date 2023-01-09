#include "RegistrationPointsVTK.hpp"
#include "FindRegistrationPoints.hpp"
#include "LeastSquaresScaleICP.hpp"
#include "TemplatePoints.hpp"
#include "ImplantTools.hpp"

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
    //data = new RegistrationPointsVTK(pKnee);
}

FindRegistrationPoints::~FindRegistrationPoints()
{
    /*data = NULL;
    delete data;*/
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
        myOfficialPoints = templateObj.mOfficialLeftPoints;
        myCheckPoints = templateObj.mLeftCheckPoints;
    }
    else
    {
        sourceAP = templateObj.vectorRightAP;
        sourceTEA = templateObj.vectorRightTEA;
        sourceCenter = templateObj.centerRight;
        tScale = (targetSize / templateObj.sizeLeft);
        myTemplatePoints = templateObj.mTemplateRight;
        myOfficialPoints = templateObj.mOfficialRightPoints;
        myCheckPoints = templateObj.mRightCheckPoints;
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

    /*
        std::vector<RegistrationPoints> result;

    std::vector<PointTypeITK> lateralBorder = PointToITKVector(data->GetLateralBorderPointsFemur());

    std::vector<PointTypeITK> latOut, latIn;

    for (int i = 0; i < lateralBorder.size(); i++)
    {
        if (i % 2 == 0)
        {
            latOut.push_back(lateralBorder[i]);
        }
        else
        {
            latIn.push_back(lateralBorder[i]);
        }
    }

    std::vector<PointTypeITK> medialBorder = PointToITKVector(data->GetMedialBorderPointsFemur());

    std::vector<PointTypeITK> medOut, medIn;

    for (int i = 0; i < medialBorder.size(); i++)
    {
        if (i % 2 == 0)
        {
            medOut.push_back(medialBorder[i]);
        }
        else
        {
            medIn.push_back(medialBorder[i]);
        }
    }

    std::vector<PointTypeITK> lateralOblique = PointToITKVector(data->GetLateralObliquePointsFemur());

    std::vector<PointTypeITK> sagital = PointToITKVector(data->GetSagitalPointsFemur());

    result.push_back(RegistrationPoints(latOut));

    result.push_back(RegistrationPoints(latIn));

    result.push_back(RegistrationPoints(lateralOblique));

    result.push_back(RegistrationPoints(sagital));

    result.push_back(RegistrationPoints(medIn));

    result.push_back(RegistrationPoints(medOut));

    return result;
    */
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
    std::vector<Point> myInterior, myExterior;
    std::vector<Point> myCheckPoints;

    if (mKnee.getIsRight() == false)
    {
        sourceAP = templateObj.vectorLeftAP;
        sourceTEA = templateObj.vectorLeftTEA;
        sourceCenter = templateObj.centerLeft;
        tScale = ((axisSize.first / templateObj.sizeLeftAP) + (axisSize.second / templateObj.sizeLeftTEA)) / 2.0;
        myTemplatePoints = templateObj.mTemplateLeft;
        myOfficialPoints = templateObj.mOfficialLeftPoints;

        myInterior = templateObj.mOfficialLeftPointsInterior;
        myExterior = templateObj.mOfficialLeftPointsExterior;
        myCheckPoints = templateObj.mLeftCheckPoints;
    }
    else
    {
        sourceAP = templateObj.vectorRightAP;
        sourceTEA = templateObj.vectorRightTEA;
        sourceCenter = templateObj.centerRight;
        tScale = ((axisSize.first / templateObj.sizeRightAp) + (axisSize.second / templateObj.sizeRightTea)) / 2.0;
        myTemplatePoints = templateObj.mTemplateRight;
        myOfficialPoints = templateObj.mOfficialRightPoints;

        myInterior = templateObj.mOfficialRightPointsInterior;
        myExterior = templateObj.mOfficialRightPointsExterior;
        myCheckPoints = templateObj.mRightCheckPoints;
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

    for (int i = 0; i < myCheckPoints.size(); i++)
    {
        cv::Mat transformPointMat = scale * (myRotation * (myCheckPoints[i]).ToMatPoint()) + myTranslation;
        Point transformPoint = Point(transformPointMat);
        double myClosest[3];
        double pnt[3] = { transformPoint.x, transformPoint.y, transformPoint.z };
        implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);
        pCheckPoints.push_back(Point(myClosest[0], myClosest[1], myClosest[2]));
    }

    /*std::vector<Point> lateralUp, medialUp, transversePointsMedial;

    data->getMedialPlateauCirclePoints(medialUp, transversePointsMedial);
    data->getLateralPlateauCirclePoints(lateralUp);

    result.push_back(RegistrationPoints(lateralUp));
    result.push_back(RegistrationPoints(medialUp));
    result.push_back(RegistrationPoints(transversePointsMedial));

    result.push_back(RegistrationPoints(data->GetTibiaLateralObliquePoints()));
    result.push_back(RegistrationPoints(data->GetTibiaMedialObliquePoints()));
    result.push_back(RegistrationPoints(data->GetTibiaOnTuberAxisPoints()));
    result.push_back(RegistrationPoints(data->GetTibiaOnTuberSidePoints()));*/

    return result;
}

/*
double FindRegistrationPoints::RegisterFemurTest(const std::vector<cv::Point3d>& templatePoints, const Point& vectorTea, const Point& vectorAp, const Point& center, double pScale)
{
    Point targetAP, targetTEA, targetAxis, targetCenter;
    targetAP = mKnee.getFemurDirectVectorAP();
    targetTEA = mKnee.getFemurVectorLateralTEA();
    targetAxis = targetAP.cross(targetTEA);
    targetAxis.normalice();
    targetCenter = (mKnee.getLateralEpicondyle() + mKnee.getMedialEpicondyle()) / 2.0;

    Point sourceAP, sourceTEA, sourceAxis, sourceCenter;
    sourceAP = vectorAp;
    sourceTEA = vectorTea;
    sourceAxis = sourceAP.cross(sourceTEA);
    sourceCenter = pScale * center;

    sourceAP.normalice();
    sourceTEA.normalice();
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

    data.at<double>(6, 0) = 1.0;

    std::vector<cv::Point3d> myTemplatePoints;
    for (int i = 0; i < templatePoints.size(); i++)
    {
        cv::Point3d temp = pScale * templatePoints[i];
        myTemplatePoints.push_back(temp);
    }

    LeastSquaresScaleICP registerObj(myTemplatePoints);

    double result = registerObj.LeastSquaresScale(mKnee.GetFemurPoly(), data);
    std::cout << "Error final: " << result << "  scale: " << data.at<double>(6, 0) << std::endl;
    return result;
}

double FindRegistrationPoints::RegisterTibiaLeft(const std::vector<cv::Point3d>& templatePoints, const Point& vectorTea, const Point& vectorAp, const Point& center, double sizeAp, double sizeTea)
{
    Point targetAP, targetTEA, targetAxis, targetCenter;
    std::pair<double, double> axisSize = mKnee.getTibiaAutomaticAxis(targetAP, targetTEA, targetCenter);
    targetAxis = targetAP.cross(targetTEA);
    targetAxis.normalice();

    double tScale = ((axisSize.first / sizeAp) + (axisSize.second / sizeTea)) / 2.0;

    Point sourceAP, sourceTEA, sourceAxis, sourceCenter;
    sourceAP = vectorAp;
    sourceTEA = vectorTea;
    sourceAxis = sourceAP.cross(sourceTEA);
    sourceCenter = center;

    sourceAP.normalice();
    sourceTEA.normalice();
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

    LeastSquaresScaleICP registerObj(templatePoints);

    double result = registerObj.LeastSquaresScale(mKnee.GetTibiaPoly(), data);
    std::cout << "Error final: " << result << "  Init scale: " << tScale << std::endl;
    return result;
}

double FindRegistrationPoints::RegisterTibiaPoints()
{
    TemplateTibia templateObj;

    Point targetAP, targetTEA, targetAxis, targetCenter;
    Point sourceAP, sourceTEA, sourceAxis, sourceCenter;
    double tScale;
    std::vector<cv::Point3d> myTemplatePoints;

    std::pair<double, double> axisSize = mKnee.getTibiaAutomaticAxis(targetAP, targetTEA, targetCenter);
    targetAxis = targetAP.cross(targetTEA);
    targetAxis.normalice();

    if (mKnee.getIsRight() == false)
    {
        sourceAP = templateObj.vectorLeftAP;
        sourceTEA = templateObj.vectorLeftTEA;
        sourceCenter = templateObj.centerLeft;
        tScale = ((axisSize.first / templateObj.sizeLeftAP) + (axisSize.second / templateObj.sizeLeftTEA)) / 2.0;
        myTemplatePoints = templateObj.mTemplateLeft;
    }
    else
    {
        sourceAP = templateObj.vectorRightAP;
        sourceTEA = templateObj.vectorRightTEA;
        sourceCenter = templateObj.centerRight;
        tScale = ((axisSize.first / templateObj.sizeRightAp) + (axisSize.second / templateObj.sizeRightTea)) / 2.0;
        myTemplatePoints = templateObj.mTemplateRight;
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

    double result = registerObj.LeastSquaresScale(mKnee.GetTibiaPoly(), data);
    std::cout << "Error final: " << result << std::endl;
    return result;
}

double FindRegistrationPoints::RegisterTibiaRight(const std::vector<cv::Point3d>& templatePoints, const Point& vectorTea, const Point& vectorAp, const Point& center, double sizeAp, double sizeTea)
{
    Point targetAP, targetTEA, targetAxis, targetCenter;
    std::pair<double, double> axisSize = mKnee.getTibiaAutomaticAxis(targetAP, targetTEA, targetCenter);
    targetAxis = targetAP.cross(targetTEA);
    targetAxis.normalice();

    double tScale = ((axisSize.first / sizeAp) + (axisSize.second / sizeTea)) / 2.0;

    Point sourceAP, sourceTEA, sourceAxis, sourceCenter;
    sourceAP = vectorAp;
    sourceTEA = vectorTea;
    sourceAxis = sourceAP.cross(sourceTEA);
    sourceCenter = center;

    sourceAP.normalice();
    sourceTEA.normalice();
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

    LeastSquaresScaleICP registerObj(templatePoints);

    double result = registerObj.LeastSquaresScale(mKnee.GetTibiaPoly(), data);
    std::cout << "Error final: " << result << "  Init scale: " << tScale << std::endl;
    return result;
}



std::vector<RegistrationPoints> FindRegistrationPoints::GetRegistrationPointsTibiaTemplate(bool isLeft)
{
    std::vector<RegistrationPoints> result;
    std::vector<Point> points;

    if (isLeft == true)
    {
        points = data->GetTibiaTemplatePoints();
    }
    else
    {
        points = data->GetTibiaTemplatePointsRight();
    }

    result.push_back(RegistrationPoints(points));

    return result;
}

std::vector<cv::Point3d> FindRegistrationPoints::GetRegistrationPointsTibiaTemplateLikeCV(bool isLeft)
{
    std::vector<Point> points;
    if (isLeft == true)
    {
        points = data->GetTibiaTemplatePoints();
    }
    else
    {
        points = data->GetTibiaTemplatePointsRight();
    }
    std::vector<cv::Point3d> result;
    for (int i = 0; i < points.size(); i++)
    {
        result.push_back(points[i]);
    }
    return result;
}

std::vector<cv::Point3d> FindRegistrationPoints::GetRegistrationPointsFemurTemplateLikeCV(bool isLeft)
{
    std::vector<Point> points;
    if (isLeft == true)
    {
        points = data->GetFemurTemplatePoints(isLeft);
    }
    else
    {
        points = data->GetFemurTemplatePoints(isLeft);
    }
    std::vector<cv::Point3d> result;
    for (int i = 0; i < points.size(); i++)
    {
        result.push_back(points[i]);
    }
    return result;
}
*/