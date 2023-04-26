#include "HipPoints.hpp"


HipPoints::HipPoints(const vtkSmartPointer<vtkPolyData>& pHip, const Point& pTopPoint, const Plane& pObliqueTransverse)
{
    mHip = pHip;
    mObliqueTransverse.init(pObliqueTransverse.getNormalVector(), pObliqueTransverse.getPoint());
    mObliqueTransverse.reverseByPoint(pTopPoint);
    mTop = pTopPoint;
}

void HipPoints::GetPointOnTop(const Point& a, const Point& c, std::vector<cv::Point3d>& pResult)
{
    Plane cutPlane;
    cutPlane.init(a, mTop, c);
    auto contour = ImplantTools::getMaxContour(mHip, cutPlane.getNormalVector(), cutPlane.getPoint());
    std::vector<Plane> Condition = { mObliqueTransverse };
    std::vector<Point> myResult;
    ImplantTools::GetPointsOnContourSort(contour, Condition, mObliqueTransverse, myResult);

    reduceSortPointsByRadius(myResult, 5, pResult);

    //for (int i = 0; i < myResult.size(); i++)
    //{
    //    pResult.push_back(myResult[i]);
    //}
}

void HipPoints::GetTransversal(const Point& a, const Point& b, const Point& c, std::vector<cv::Point3d>& pResult)
{
    Plane cutPlane;
    cutPlane.init(a, b, c);
    auto contour = ImplantTools::getMaxContour(mHip, cutPlane.getNormalVector(), cutPlane.getPoint());

    std::vector<Point> allPoints;

    std::list<std::pair<vtkIdType, vtkIdType>> lines;
    ImplantTools::ExtractSortLines(contour, lines);

    std::list<std::pair<vtkIdType, vtkIdType>>::iterator it11, it22;
    it11 = lines.begin();
    it22 = lines.end();

    for (; it11 != it22; ++it11)
    {
        double pnt[3];
        contour->GetPoint(it11->first, pnt);
        allPoints.push_back(Point(pnt[0], pnt[1], pnt[2]));
    }

    reduceSortPointsByRadius(allPoints, 5, pResult);
}

void HipPoints::GetSagitalDown(const Point& a, const Point& b, const Point& c, std::vector<cv::Point3d>& pResult)
{
    Plane cutPlane;
    cutPlane.init(a, b, c);
    auto contour = ImplantTools::getMaxContour(mHip, cutPlane.getNormalVector(), cutPlane.getPoint());

    std::list<std::pair<vtkIdType, vtkIdType>> lines;
    ImplantTools::ExtractSortLines(contour, lines);

    std::list<std::pair<vtkIdType, vtkIdType>>::iterator it1, it2;
    it1 = lines.begin();
    it2 = lines.end();

    vtkIdType nearA, nearC;

    nearA = ImplantTools::GetNearestPoints(contour, a);
    nearC = ImplantTools::GetNearestPoints(contour, c);

    std::vector<vtkIdType> allPoints;
    int cont = 0;
    int pos = -1;

    for (; it1 != it2; ++it1)
    {
        allPoints.push_back(it1->first);
        if (it1->first == nearA)
        {
            pos = cont;
        }
        cont++;
    }

    int tSize = cont;

    if (lines.back().second == nearA && pos == -1)
    {
        pos = cont - 1;
    }

    if (pos > 0)
    {
        std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());
    }
    else if (pos == -1)
    {
        return;
    }

    double pnt[3];
    contour->GetPoint(allPoints[5], pnt);
    Point nearPoint(pnt[0], pnt[1], pnt[2]);

    contour->GetPoint(allPoints[allPoints.size() - 5], pnt);
    Point farPoint(pnt[0], pnt[1], pnt[2]);

    if (ImplantTools::getDistanceBetweenPoints(nearPoint, c) > ImplantTools::getDistanceBetweenPoints(farPoint, c))
    {
        std::reverse(allPoints.begin(), allPoints.end());
    }

    std::vector<Point> myResult;

    for (int i = 0; i < allPoints.size(); i++)
    {
        double pnt[3];
        contour->GetPoint(allPoints[i], pnt);
        myResult.push_back(Point(pnt[0], pnt[1], pnt[2]));
        if (allPoints[i] == nearC)
        {
            break;
        }
    }

    reduceSortPointsByRadius(myResult, 4, pResult);
}

void HipPoints::reduceSortPointsByRadius(const std::vector<Point>& sortPoints, double radius, std::vector<cv::Point3d>& result)
{
    int tSize = sortPoints.size();

    if (tSize < 3)
    {
        return;
    }

    result.push_back(sortPoints[0]);

    for (int i = 1; i < tSize; i++)
    {
        if (ImplantTools::isPointInsideSphere(result[result.size() - 1], radius, sortPoints[i]) == false)
        {
            Line myLine = Line::makeLineWithPoints(sortPoints[i], sortPoints[i - 1]);
            std::pair<Point, Point> intercep;

            myLine.getInterceptionSphere(result[result.size() - 1], radius, intercep);

            if (ImplantTools::getDistanceBetweenPoints(sortPoints[i], intercep.first, true) < ImplantTools::getDistanceBetweenPoints(sortPoints[i], intercep.second, true))
            {
                result.push_back(intercep.first);
            }
            else
            {
                result.push_back(intercep.second);
            }
        }

    }
}