#include "PatellaImplantMatch.hpp"
#include "ImplantsException.hpp"
#include "vtkCutter.h"
#include "vtkPlane.h"
#include "ImplantTools.hpp"
#include "vtkPlaneCollection.h"
#include "vtkClipClosedSurface.h"
#include "ConvexHull.hpp"


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


PatellaImplantMatch::PatellaImplantMatch()
{
    isInit = false;
}

void PatellaImplantMatch::init(const PatellaImplant& implant, const Knee& knee)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_PATELLA_IMPLANT_MATCH;
    }
    this->implant = implant;
    this->knee = knee;
    makeRotationMatrix();
    makeTranslationMatrix();
    isInit = true;
}

Plane PatellaImplantMatch::transformPlane(const Plane& plane, const itk::Rigid3DTransform<>::Pointer pTransform) const
{
    cv::Mat rotation = Rigid3DTransformToCVRotation(pTransform);
    cv::Mat translation = Rigid3DTransformToCVTranslation(pTransform);

    cv::Mat transformNormalVector = rotation * plane.getNormalVectorMat();
    cv::Mat transformPoint = (rotation * plane.getPointMat()) + translation;
    Plane transformPlane;
    transformPlane.init(Point(transformNormalVector), Point(transformPoint));
    return transformPlane;
}

void PatellaImplantMatch::makeRotationMatrix()
{
    Point fromVector = implant.getNormalVector();
    Point toVector = knee.getPatellaFrontPlane().getNormalVector();

    rotationMatrix = ImplantTools::GetGeneralRotateTransformVectors(fromVector, toVector);
}

void PatellaImplantMatch::makeTranslationMatrix()
{
    Point patellaCenter = knee.getPatellaDistalPosteriorPoint() - (implant.getThickness() * knee.getPatellaFrontPlane().getNormalVector());
    translationMatrix = patellaCenter.ToMatPoint() - (rotationMatrix * implant.getCentralPointOnBase().ToMatPoint());
}

Plane PatellaImplantMatch::transformPlane(const Plane& plane) const
{
    cv::Mat transformNormalVector = rotationMatrix * plane.getNormalVectorMat();
    cv::Mat transformPoint = (rotationMatrix * plane.getPointMat()) + translationMatrix;
    Plane transformPlane;
    transformPlane.init(Point(transformNormalVector), Point(transformPoint));
    return transformPlane;
}

Plane PatellaImplantMatch::getPatellaCutPlane() const
{
    return transformPlane(implant.getBasePlane());
}

Point PatellaImplantMatch::transformImplantPoint(const Point& pPoint) const
{
    cv::Mat mat(3, 1, CV_64FC1);
    mat.at <double>(0, 0) = pPoint.x;
    mat.at <double>(1, 0) = pPoint.y;
    mat.at <double>(2, 0) = pPoint.z;

    cv::Mat transformPoint = (rotationMatrix * mat) + translationMatrix;
    return Point(transformPoint);
}

itk::Matrix< double, 3, 3 > PatellaImplantMatch::GetRotationMatrix() const
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

itk::Vector< double, 3 > PatellaImplantMatch::GetTranslationMatrix() const
{
    itk::Vector< double, 3 > translation;
    translation[0] = translationMatrix.at <double>(0, 0);
    translation[1] = translationMatrix.at <double>(1, 0);
    translation[2] = translationMatrix.at <double>(2, 0);

    return translation;
}

std::vector<PointTypeITK> PatellaImplantMatch::GetHullPoints(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut, double distance, int amount) const
{
    std::vector<PointTypeITK> hull, tempHull;
    Plane myPlane = transformPlane(implant.getBasePlane(), pTransformIn);
    myPlane.reverseByPoint(implant.getTopPoint(), false);

    Point sagitalVector = -knee.getPatella().getPatellaInferiorVector();

    cv::Mat myRotation = getTransformToRobot(myPlane, sagitalVector);

    vtkSmartPointer<vtkPolyData> contourMax = ImplantTools::getContours(knee.GetPatellaPoly(), myPlane.getNormalVector(), myPlane.getPoint());

    vtkSmartPointer<vtkPoints> vtkMyPoints = contourMax->GetPoints();
    int tVtkPointSize = vtkMyPoints->GetNumberOfPoints();

    std::vector<Point> contourPoints;

    for (int i = 0; i < tVtkPointSize; i++)
    {
        double pnt[3];
        contourMax->GetPoint(i, pnt);
        Point currentPoint(pnt[0], pnt[1], pnt[2]);
        contourPoints.push_back(currentPoint);
    }

    if (contourPoints.size() < 20)
    {
        throw ImplantExceptionCode::NOT_ENOUGH_POINTS_ON_CUT_PLANE;
    }

    ConvexHull hullObj(contourPoints, myRotation);
    std::vector<Point> hullPoints = hullObj.GetConvexHull();

    if (hullPoints.size() < 3)
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_CONVEX_HULL_ON_PATELLA_CUT_PLANE;
    }

    std::vector<Point> concaveHull;

    Point computePoint, midPoint;

    for (int i = 0; i < hullPoints.size() - 1; i++)
    {
        midPoint = midPoint + hullPoints[i];

        computePoint = movePointAtNormal(hullPoints[i], hullPoints[i + 1], myRotation, distance);
        concaveHull.push_back(computePoint);

        if (i == hullPoints.size() - 2)
        {
            midPoint = midPoint + hullPoints[i + 1];

            computePoint = movePointAtNormal(hullPoints[i], hullPoints[i + 1], myRotation, distance, true);
            concaveHull.push_back(computePoint);
        }
    }

    midPoint = midPoint / double(hullPoints.size());

    tempHull = increaseVectorToAmount(concaveHull, amount);

    double angle = (40.0) * PI / 180.0;

    Point latVector = myPlane.getProjectionVector(knee.getPatella().getPatellaLateralVector());
    Point infVector = myPlane.getProjectionVector(knee.getPatella().getPatellaInferiorVector());
    latVector.normalice();
    infVector.normalice();

    Point pointLat = midPoint + latVector + tan(angle) * infVector;
    Point pointMed = midPoint - latVector + tan(angle) * infVector;

    Plane obliqueLatPlaneDown = myPlane.getPerpendicularPlane(midPoint, pointLat);
    Plane obliqueMedPlaneDown = myPlane.getPerpendicularPlane(midPoint, pointMed);

    obliqueLatPlaneDown.reverseByPoint(pointMed);
    obliqueMedPlaneDown.reverseByPoint(pointLat);

    int initPos = -1;

    for (int i = 0; i < tempHull.size(); i++)
    {
        Point temp = Point(tempHull[i][0], tempHull[i][1], tempHull[i][2]);

        if (obliqueLatPlaneDown.eval(temp) >= 0 && obliqueMedPlaneDown.eval(temp) >= 0)
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
        std::rotate(tempHull.begin(), tempHull.begin() + initPos, tempHull.end());
    }

    std::vector<Point> finalPoints;

    for (int i = 0; i < tempHull.size(); i++)
    {
        Point temp = Point(tempHull[i][0], tempHull[i][1], tempHull[i][2]);

        if (!(obliqueLatPlaneDown.eval(temp) >= 0 && obliqueMedPlaneDown.eval(temp) >= 0))
        {
            finalPoints.push_back(temp);
        }
    }

    hull = increaseVectorToAmount(finalPoints, amount);

    itk::Vector< double, 3 > translate;

    Point tTemp = myRotation * (midPoint.ToMatPoint());
    translate[0] = -tTemp.x;
    translate[1] = -tTemp.y;
    translate[2] = -tTemp.z;

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

Point PatellaImplantMatch::movePointAtNormal(const Point& movePoint, const Point& nextPoint, const cv::Mat& rotationZ, double distance, bool changeMove, bool clockWise) const
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

std::vector<PointTypeITK> PatellaImplantMatch::increaseVectorToAmount(const std::vector<Point>& points, int amount) const
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

vtkSmartPointer<vtkPolyData> PatellaImplantMatch::GetCuttingPatella() const
{
    Plane PlaneA = transformPlane(implant.getBasePlane());
    //Point ref = knee.getFemurKneeCenter() + (10000.0 * knee.getDirectVectorFemurAxis()) - (10000.0 * knee.getFemurDirectVectorAP());
    PlaneA.reverse();

    cv::Point3d planeNormal, planePoint;

    vtkNew<vtkPlane> vtkPlane;

    planeNormal = PlaneA.getNormalVector();
    planePoint = PlaneA.getPoint();
    vtkPlane->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
    vtkPlane->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);

    vtkNew<vtkPlaneCollection> patellaPlanes;
    patellaPlanes->AddItem(vtkPlane);

    vtkNew<vtkClipClosedSurface> patellaClipper;
    patellaClipper->SetInputData(knee.GetPatellaPoly());
    patellaClipper->SetClippingPlanes(patellaPlanes);
    patellaClipper->Update();

    return patellaClipper->GetOutput();
}