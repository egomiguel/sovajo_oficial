#include "ImplantsException.hpp"
#include "Patella.hpp"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkCleanPolyData.h"
#include "ImplantTools.hpp"

using namespace TKA::IMPLANTS;

Patella::Patella()
{
    isInit = false;
}

void Patella::init(const Point& pRotationPoint, const Point& pLateralPoint, const Point& pMedialPoint, const Point& pInferiorPoint, const vtkSmartPointer<vtkPolyData> pPatellaPoly)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_PATELLA;
    }

    vtkNew<vtkCleanPolyData> CleanPatella;

    vtkNew<vtkPolyDataConnectivityFilter> patellaConnectivityFilter;
    patellaConnectivityFilter->SetInputData(pPatellaPoly);
    patellaConnectivityFilter->SetExtractionModeToLargestRegion();
    patellaConnectivityFilter->Update();

    CleanPatella->SetInputData(patellaConnectivityFilter->GetOutput());
    CleanPatella->Update();

    mRotationPoint = pRotationPoint;
    mLateralPoint = pLateralPoint;
    mMedialPoint = pMedialPoint;
    mInferiorPoint = pInferiorPoint;
    mPatellaPoly = CleanPatella->GetOutput();

    isInit = true;
}

bool Patella::getIsInit()
{
    return isInit;
}

Plane Patella::getPatellaPlane(bool isRight, Point& furthestPoint) const
{
    if (isInit == false)
    {
        return Plane();
    }

    Plane planeFront = getFrontPlane(isRight);

    if (planeFront.getIsInit() == false)
    {
        throw ImplantExceptionCode::CAN_NOT_FIT_PATELLA_FRONTAL_PLANE;
    }

    double distance = 0, temp = 0;
    vtkSmartPointer<vtkPoints> points = mPatellaPoly->GetPoints();
    int tSize = points->GetNumberOfPoints();
    Point farPoint;
    bool isOk = false;
    for (int i = 0; i < tSize; i++)
    {
        double pnt[3];
        points->GetPoint(i, pnt);
        temp = planeFront.eval(pnt);

        if (temp > distance)
        {
            distance = temp;
            farPoint.x = pnt[0];
            farPoint.y = pnt[1];
            farPoint.z = pnt[2];
            isOk = true;
        }
    }

    if (isOk == false)
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_PATELLA_DISTAL_POSTERIOR_POINT;
    }

    Plane tempPlane;
    tempPlane.init(planeFront.getNormalVector(), farPoint);

    furthestPoint = tempPlane.getProjectionPoint(getPatellaCenter());

    return planeFront;
}

Plane Patella::getFrontPlane(bool isRight) const
{
    Plane refPlane;
    refPlane.init(mLateralPoint, mMedialPoint, mInferiorPoint);

    Point vectorHorizontal = mLateralPoint - mMedialPoint;
    Line lineML(vectorHorizontal, mLateralPoint);
    Point supPoint = lineML.getProjectPoint(mInferiorPoint);

    Point vectorSupInf = mInferiorPoint - supPoint;
    Point vectorRef;

    if (isRight == true)
    {
        vectorRef = mMedialPoint - mLateralPoint;

    }
    else
    {
        vectorRef = mLateralPoint - mMedialPoint;
    }

    Point vectorGood = vectorRef.cross(vectorSupInf);
    vectorGood.normalice();
    refPlane.reverseByNormal(vectorGood);

    double box[6];
    mPatellaPoly->GetBounds(box);

    std::vector<double> dimesions = { abs(box[1] - box[0]), abs(box[3] - box[2]), abs(box[5] - box[4]) };
    std::sort(dimesions.begin(), dimesions.end());

    refPlane.movePlaneOnNormal(-dimesions[0]);

    Point center = refPlane.getProjectionPoint(getPatellaCenter());
    Point pointLat = refPlane.getProjectionPoint(mLateralPoint);
    Point pointMed = refPlane.getProjectionPoint(mMedialPoint);
    Point pointInf = refPlane.getProjectionPoint(mInferiorPoint);

    pointLat = pointLat + 0.1 * (center - pointLat);
    pointMed = pointMed + 0.1 * (center - pointMed);
    pointInf = pointInf + 0.1 * (center - pointInf);

    std::vector<Point> result;

    vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
    implicitPolyDataDistance->SetInput(mPatellaPoly);

    getPatellaFront(pointLat, pointInf, center, implicitPolyDataDistance, result);
    getPatellaFront(pointMed, pointInf, center, implicitPolyDataDistance, result);

    bool isOk;
    Plane planeFront = Plane::getBestPlane(result, isOk);

    if (isOk == true)
    {
        center = getPatellaCenter() + dimesions[2] * refPlane.getNormalVector();
        planeFront.reverseByPoint(center);

        double angle = ImplantTools::getAngleBetweenVectorsDegree(planeFront.getNormalVector(), refPlane.getNormalVector());

        if (angle > 25)
        {
            planeFront.deletePlane();
            planeFront.init(refPlane.getNormalVector(), refPlane.getPoint());
        }
    }

    return planeFront;
}

Point Patella::getPatellaCenter() const
{
    if (isInit == false)
    {
        return Point();
    }

    double pnt[3];
    mPatellaPoly->GetCenter(pnt);
    return Point(pnt[0], pnt[1], pnt[2]);
}

double Patella::getPatellaDiameter() const
{
    if (isInit == false)
    {
        return -1;
    }

    double box[6];
    mPatellaPoly->GetBounds(box);

    std::vector<double> dimesions = { abs(box[1] - box[0]), abs(box[3] - box[2]), abs(box[5] - box[4]) };
    std::sort(dimesions.begin(), dimesions.end());
    double middleSize = dimesions[1];
    return middleSize;
}

double Patella::getPatellaThickness() const
{
    if (isInit == false)
    {
        return -1;
    }

    double box[6];
    mPatellaPoly->GetBounds(box);

    std::vector<double> dimesions = { abs(box[1] - box[0]), abs(box[3] - box[2]), abs(box[5] - box[4]) };
    std::sort(dimesions.begin(), dimesions.end());
    double thickness = dimesions[0];
    return thickness;
}

Point Patella::getPatellaInferiorVector() const
{
    if (isInit == false)
    {
        return Point();
    }

    Plane refPlane;
    refPlane.init(mLateralPoint, mMedialPoint, mInferiorPoint);
    Point vector = refPlane.getProjectionPoint(mInferiorPoint) - refPlane.getProjectionPoint(getPatellaCenter());
    vector.normalice();
    return vector;
}

Point Patella::getPatellaLateralVector() const
{
    if (isInit == false)
    {
        return Point();
    }

    Plane refPlane;
    refPlane.init(mLateralPoint, mMedialPoint, mInferiorPoint);

    Line refLine = Line::makeLineWithPoints(refPlane.getProjectionPoint(mInferiorPoint), refPlane.getProjectionPoint(getPatellaCenter()));
    Point projLat = refPlane.getProjectionPoint(mLateralPoint);
    Point projLine = refLine.getProjectPoint(projLat);

    Point vector = projLat - projLine;
    vector.normalice();
    return vector;
}

Point Patella::getPatellaRotationPoint() const
{
    if (isInit == false)
    {
        return Point();
    }

    return mRotationPoint;
}

vtkSmartPointer<vtkPolyData> Patella::getPatellaPoly() const
{
    if (isInit == false)
    {
        return nullptr;
    }

    return mPatellaPoly;
}

void Patella::getPatellaFront(const Point& a, const Point& b, const Point& center, const vtkSmartPointer<vtkImplicitPolyDataDistance>& polyDistance, std::vector<Point>& points) const
{
    std::vector<Point> vectorTemp;
    double coef = 0;

    for (double i = 0; i <= 10; i++)
    {
        coef = i * (0.1);
        Point temp = a + coef * (b - a);
        vectorTemp.push_back(temp);
    }

    Line myLine = Line::makeLineWithPoints(a, b);
    double distance = myLine.getDistanceFromPoint(center);
    Point direction = center - (myLine.getProjectPoint(center));
    direction.normalice();

    coef = distance / 10.;

    for (int i = 0; i < vectorTemp.size(); i++)
    {
        for (double j = 0; j < 10; j++)
        {
            Point temp = vectorTemp[i] + j * coef * direction;

            double myClosest[3];

            double pnt[3] = { temp.x, temp.y, temp.z };
            polyDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);

            points.push_back(Point(myClosest[0], myClosest[1], myClosest[2]));
        }
    }
}




