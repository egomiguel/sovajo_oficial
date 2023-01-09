#include "GenerateSurfaceOld.hpp"
#include <list>
#include "vtkTriangle.h"
#include "vtkKdTreePointLocator.h"
#include "SegmentationException.hpp"
#include "vtkPolyLine.h"
#include "vtkCellLocator.h"
#include "vtkRuledSurfaceFilter.h"
#include "vtkContourTriangulator.h"
#include "vtkAppendPolyData.h"

const double EPSILON = std::numeric_limits<double>::epsilon();

const int lessAmount = 10;

const double closeAngle = 100.0 * vtkMath::Pi() / 180.0;

vtkIdType FixPos(vtkIdType pos, vtkIdType size, vtkIdType offset, vtkIdType rot)
{
    if (rot == 0)
    {
        return pos;
    }

    vtkIdType temp = pos + rot - offset;
    if (temp < size)
    {
        return temp + offset;
    }
    else
    {
        return temp - size + offset;
    }
}

SPlane GenerateSurfaceOld::GetFitPlane(const vtkSmartPointer<vtkPoints> pointsList) const
{
    double p1[3], p2[3];
    vtkIdType size = pointsList->GetNumberOfPoints();

    pointsList->GetPoint(0, p1);
    pointsList->GetPoint(size / 3, p2);

    std::vector<cv::Point3d> points;

    for (int i = (size / 3) + 1; i < size; i++)
    {
        double p3[3];
        pointsList->GetPoint(i, p3);
        if (AreCollinear(p1, p2, p3) == false)
        {
            points.push_back(VtkToCV(p3));
        }
    }

    cv::Point3d cvP1, cvP2, directVector1, directVector2;
    cv::Point3d normal = cv::Point3d(0, 0, 0);
    cvP1 = VtkToCV(p1);
    cvP2 = VtkToCV(p2);

    for (int i = 0; i < points.size(); i++)
    {
        directVector1 = cvP1 - points[i];
        directVector2 = cvP2 - points[i];
        normal = normal + (directVector1.cross(directVector2));
    }
    double lenght = double(points.size());
    normal = normal / lenght;
    SPlane plane;
    plane.init(normal, cvP1);
    return plane;
}

GenerateSurfaceOld::GenerateSurfaceOld(std::vector<vtkSmartPointer<vtkPolyData>> pSlices)
{
    slices = pSlices;
    allPoints = vtkSmartPointer<vtkPoints>::New();
    cells = vtkSmartPointer<vtkCellArray>::New();
}

double GenerateSurfaceOld::GetDistance(const double P1[3], const double P2[3], bool square) const
{
    double a, b, c;
    a = P1[0] - P2[0];
    b = P1[1] - P2[1];
    c = P1[2] - P2[2];

    if (square == true)
    {
        return (a*a + b * b + c * c);
    }
    else
    {
        return sqrt(a*a + b * b + c * c);
    }

}

double GenerateSurfaceOld::GetDistance(const cv::Point3d& P1, const cv::Point3d& P2) const
{
    cv::Point3d diff = P2 - P1;
    return sqrt(diff.dot(diff));
}

double GenerateSurfaceOld::ExtractSortLines(const vtkSmartPointer<vtkPolyData> polyData, std::list<std::pair<vtkIdType, vtkIdType>>& lines) const
{
    vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
    vtkSmartPointer<vtkCellArray> cells = polyData->GetLines();
    cells->InitTraversal();
    vtkNew<vtkIdList> myList;
    std::list<std::pair<vtkIdType, vtkIdType>> restLines;
    double lessDistance = 9999999.0;
    double totalDistance = 0;
    while (cells->GetNextCell(myList))
    {
        vtkIdType id1 = myList->GetId(0);
        vtkIdType id2 = myList->GetId(1);

        double p1[3], p2[3];
        points->GetPoint(id1, p1);
        points->GetPoint(id2, p2);
        double tempDistance = GetDistance(p1, p2);
        totalDistance = totalDistance + tempDistance;

        if (tempDistance < lessDistance)
        {
            lessDistance = tempDistance;
        }

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

            double p1[3], p2[3];
            points->GetPoint(id1, p1);
            points->GetPoint(id2, p2);
            double tempDistance = GetDistance(p1, p2);
            totalDistance = totalDistance + tempDistance;

            if (tempDistance < lessDistance)
            {
                lessDistance = tempDistance;
            }

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
    return totalDistance / double(lines.size());
}

double GenerateSurfaceOld::ExtractPoints(const vtkSmartPointer<vtkPolyData> polyData, vtkSmartPointer<vtkPoints> pPoints) const
{
    vtkSmartPointer<vtkPoints> points = polyData->GetPoints();

    std::list<std::pair<vtkIdType, vtkIdType>> lines;

    std::list<cv::Point3d> sortPoints, linesOut;

    double lessDistance = ExtractSortLines(polyData, lines);
    if (lessDistance > 0.5)
    {
        lessDistance = 0.5;
    }

    auto it1 = lines.begin();
    auto it2 = lines.end();

    for (; it1 != it2; ++it1)
    {
        double p[3];
        points->GetPoint((*it1).first, p);
        sortPoints.push_back(cv::Point3d(p[0], p[1], p[2]));
    }
    double totalDistance = Resample(lessDistance, sortPoints, linesOut);

    auto cvIt1 = linesOut.begin();
    auto cvIt2 = linesOut.end();

    for (; cvIt1 != cvIt2; ++cvIt1)
    {
        double p[3];
        p[0] = (*cvIt1).x;
        p[1] = (*cvIt1).y;
        p[2] = (*cvIt1).z;
        pPoints->InsertNextPoint(p);
    }

    return totalDistance;
}

std::list<cv::Point3d> GenerateSurfaceOld::ExtractPointsList(const vtkSmartPointer<vtkPolyData> polyData) const
{
    vtkSmartPointer<vtkPoints> points = polyData->GetPoints();

    std::list<std::pair<vtkIdType, vtkIdType>> lines;

    double lessDistance = ExtractSortLines(polyData, lines);

    auto it1 = lines.begin();
    auto it2 = lines.end();

    std::list<cv::Point3d> result;

    for (; it1 != it2; ++it1)
    {
        double p[3];
        points->GetPoint((*it1).first, p);
        result.push_back(cv::Point3d(p[0], p[1], p[2]));
    }

    return result;
}

double GenerateSurfaceOld::Resample(double lessDistance, const std::list<cv::Point3d>& linesIn, std::list<cv::Point3d>& linesOut) const
{
    rsize_t size = linesIn.size();
    auto it1 = linesIn.begin();
    auto it2 = std::next(it1, size - 1);
    double totalDistance = 0.0;

    for (; it1 != it2; ++it1)
    {
        linesOut.push_back(*it1);

        cv::Point3d cvP1 = *it1;
        cv::Point3d cvP2 = *(std::next(it1, 1));

        double distance = GetDistance(cvP1, cvP2);
        totalDistance = totalDistance + distance;
        double quotient = distance / lessDistance;
        if (quotient > 1.8)
        {
            cv::Point3d vector = cvP2 - cvP1;
            vector = vector / (sqrt(vector.dot(vector)));
            for (double i = lessDistance; i <= distance - lessDistance; i += lessDistance)
            {
                cv::Point3d result = cvP1 + i * vector;
                linesOut.push_back(result);

            }
        }

        if (std::next(it1, 1) == it2)
        {
            linesOut.push_back(*it2);

            cvP1 = *it2;
            cvP2 = *(linesIn.begin());

            distance = GetDistance(cvP1, cvP2);
            totalDistance = totalDistance + distance;
            quotient = distance / lessDistance;
            if (quotient > 1.8)
            {
                cv::Point3d vector = cvP2 - cvP1;
                vector = vector / (sqrt(vector.dot(vector)));
                for (double i = lessDistance; i <= distance - lessDistance; i += lessDistance)
                {
                    cv::Point3d result = cvP1 + i * vector;
                    linesOut.push_back(result);

                }
            }

        }
    }
    //std::cout << "Resample size before: " << linesOut.size() << std::endl;
    if (linesOut.size() <= lessAmount)
    {
        std::list<cv::Point3d> tempList;
        size = linesOut.size();

        it1 = linesOut.begin();
        it2 = linesOut.end();

        cv::Point3d cvP1, cvP2;
        for (; it1 != it2; ++it1)
        {
            tempList.push_back(*it1);
            cvP1 = *it1;

            if (std::next(it1, 1) == it2)
            {
                cvP2 = *(linesOut.begin());
            }
            else
            {
                cvP2 = *(std::next(it1, 1));
            }

            if (size >= 7)
            {
                tempList.push_back(cvP1 + (0.5)*(cvP2 - cvP1));
            }
            else if (size >= 5)
            {
                tempList.push_back(cvP1 + (0.33)*(cvP2 - cvP1));
                tempList.push_back(cvP1 + (0.66)*(cvP2 - cvP1));
            }
            else
            {
                tempList.push_back(cvP1 + (0.25)*(cvP2 - cvP1));
                tempList.push_back(cvP1 + (0.5)*(cvP2 - cvP1));
                tempList.push_back(cvP1 + (0.75)*(cvP2 - cvP1));
            }
        }

        linesOut.clear();
        linesOut = tempList;

    }
    //std::cout << "Resample size after: " << linesOut.size() << std::endl;
    return totalDistance;
}

void GenerateSurfaceOld::CreateCells(vtkIdType p1, vtkIdType p2, vtkIdType p3)
{
    vtkNew<vtkTriangle> triangle;
    triangle->GetPointIds()->SetId(0, p1);
    triangle->GetPointIds()->SetId(1, p2);
    triangle->GetPointIds()->SetId(2, p3);
    cells->InsertNextCell(triangle);
}

vtkIdType GenerateSurfaceOld::FindClosestPoint(const vtkSmartPointer<vtkPoints> fitPoints, const vtkSmartPointer<vtkPoints> rotPoints, vtkIdType& rotFit) const
{
    vtkNew<vtkPolyData> polydata;
    polydata->SetPoints(rotPoints);

    vtkNew<vtkKdTreePointLocator> kDTree;
    kDTree->SetDataSet(polydata);
    kDTree->BuildLocator();

    vtkIdType size = fitPoints->GetNumberOfPoints();

    vtkIdType iD, tempId;
    double distance, tempDistance;
    distance = 9999999.0;

    for (vtkIdType i = 0; i < size; i++)
    {
        double pnt[3];
        fitPoints->GetPoint(i, pnt);
        tempId = kDTree->FindClosestPoint(pnt);

        tempDistance = GetDistance(rotPoints->GetPoint(tempId), pnt, true);
        if (tempDistance < distance)
        {
            distance = tempDistance;
            iD = tempId;
            rotFit = i;
        }
    }
    return iD;
}

void GenerateSurfaceOld::StitchPoints(vtkIdType bigSetEnd, int step, int rest, vtkIdType generalOffSetBig, vtkIdType generalOffSetSmoll, vtkIdType fitGeneralOffset, vtkIdType fitSize, vtkIdType fitRot, bool localChangeOffset)
{
    int cont = 0;
    if (rest == 0)
    {
        cont = 1;
    }

    vtkIdType currentSmollPos = 0;
    vtkIdType newPos1, newPos2;
    vtkIdType bigOffSet, smollOffset;
    if (localChangeOffset == true)
    {
        bigOffSet = generalOffSetSmoll;
        smollOffset = generalOffSetBig;
    }
    else
    {
        bigOffSet = generalOffSetBig;
        smollOffset = generalOffSetSmoll;
    }
    for (vtkIdType i = 0; i < bigSetEnd - 1; i++)
    {
        if (cont < step)
        {
            if (localChangeOffset == true)
            {
                newPos1 = FixPos(currentSmollPos + smollOffset, fitSize, fitGeneralOffset, fitRot);

                CreateCells(i + bigOffSet, i + 1 + bigOffSet, newPos1);
            }
            else
            {
                newPos1 = FixPos(i + bigOffSet, fitSize, fitGeneralOffset, fitRot);
                newPos2 = FixPos(i + 1 + bigOffSet, fitSize, fitGeneralOffset, fitRot);

                CreateCells(newPos1, newPos2, currentSmollPos + smollOffset);
            }
            cont++;
        }
        else
        {
            if (localChangeOffset == true)
            {
                newPos1 = FixPos(currentSmollPos + smollOffset, fitSize, fitGeneralOffset, fitRot);
                newPos2 = FixPos(currentSmollPos + 1 + smollOffset, fitSize, fitGeneralOffset, fitRot);

                CreateCells(i + bigOffSet, newPos1, newPos2);

                newPos1 = FixPos(currentSmollPos + 1 + smollOffset, fitSize, fitGeneralOffset, fitRot);

                CreateCells(i + bigOffSet, i + 1 + bigOffSet, newPos1);
            }
            else
            {
                newPos1 = FixPos(i + bigOffSet, fitSize, fitGeneralOffset, fitRot);

                CreateCells(newPos1, currentSmollPos + smollOffset, currentSmollPos + 1 + smollOffset);

                newPos1 = FixPos(i + bigOffSet, fitSize, fitGeneralOffset, fitRot);
                newPos2 = FixPos(i + 1 + bigOffSet, fitSize, fitGeneralOffset, fitRot);

                CreateCells(newPos1, newPos2, currentSmollPos + 1 + smollOffset);
            }
            currentSmollPos++;
            rest--;
            if (rest < 1)
            {
                cont = 1;
            }
            else
            {
                cont = 0;
            }
        }
    }
    //CreateCells(totalBigSize - 1 + bigOffSet, currentSmollPos + smollOffset, initSmoll + smollOffset);
    //CreateCells(totalBigSize - 1 + bigOffSet, initBig + bigOffSet, initSmoll + smollOffset);
}

vtkSmartPointer<vtkPoints> GenerateSurfaceOld::PrepareStitchOneSlices(const vtkSmartPointer<vtkPoints> fitPoints, const vtkSmartPointer<vtkPoints> rotatePoints, bool remove)
{
    vtkIdType size1 = fitPoints->GetNumberOfPoints();
    vtkIdType size2 = rotatePoints->GetNumberOfPoints();
    vtkIdType generalOffset = allPoints->GetNumberOfPoints() - size1;
    vtkIdType offSet1, offSet2;

    if (size1 == 0 || size2 == 0)
    {
        throw SegmentationException("There are empty slices. There can be no empty slices.");
    }

    if (generalOffset < 0)
    {
        throw SegmentationException("Unable to perform slice reconstruction, check the generated points of your polydata.");
    }
    vtkIdType rotFit;
    vtkIdType nearId = FindClosestPoint(fitPoints, rotatePoints, rotFit);

    vtkNew<vtkPoints> tempRotate, tempFit;

    for (vtkIdType i = nearId; i < size2; i++)
    {
        tempRotate->InsertNextPoint(rotatePoints->GetPoint(i));
    }

    for (vtkIdType i = 0; i < nearId; i++)
    {
        tempRotate->InsertNextPoint(rotatePoints->GetPoint(i));
    }

    for (vtkIdType i = rotFit; i < size1; i++)
    {
        tempFit->InsertNextPoint(fitPoints->GetPoint(i));
    }

    for (vtkIdType i = 0; i < rotFit; i++)
    {
        tempFit->InsertNextPoint(fitPoints->GetPoint(i));
    }

    int step1, stepSize1, stepSize2, fit, rot;
    int initFit = 0;
    int initRot = 0;
    vtkIdType beginFit, endFit, beginRot, endRot;

    std::vector<int> fitSteps, rotateSteps;

    FindBestStepCombination(tempFit, tempRotate, fitSteps, rotateSteps);

    offSet1 = generalOffset;
    offSet2 = generalOffset + size1;
    for (int k = 0; k < fitSteps.size(); k++)
    {
        stepSize1 = fitSteps[k];
        stepSize2 = rotateSteps[k];
        
        fit = stepSize1 - initFit;
        rot = stepSize2 - initRot;

        MakeStitchOneSlices(fit, rot, offSet1, offSet2, generalOffset, size1, rotFit);

        initFit = stepSize1 - 1;
        initRot = stepSize2 - 1;
        offSet1 = generalOffset + stepSize1 - 1;
        offSet2 = generalOffset + size1 + stepSize2 - 1;

    }

    for (vtkIdType i = 0; i < size2; i++)
    {
        allPoints->InsertNextPoint(tempRotate->GetPoint(i));
    }

    beginFit = generalOffset + rotFit;
    if (rotFit == 0)
    {
        endFit = generalOffset + size1 - 1;
    }
    else
    {
        endFit = generalOffset + rotFit - 1;
    }

    beginRot = generalOffset + size1;
    endRot = generalOffset + size1 + size2 - 1;

    CreateCells(beginFit, endFit, beginRot);
    CreateCells(endFit, endRot, beginRot);

    return tempRotate;
}

void GenerateSurfaceOld::MakeStitchOneSlices(vtkIdType size1, vtkIdType size2, vtkIdType offSet1, vtkIdType offSet2, vtkIdType fitGeneralOffset, vtkIdType fitTotalSize, vtkIdType rotationFit)
{
    int step, rest;

    if (size1 > size2)
    {
        step = size1 / size2;
        rest = size1 % size2;
        StitchPoints(size1, step, rest, offSet1, offSet2, fitGeneralOffset, fitTotalSize, rotationFit);

    }
    else
    {
        step = size2 / size1;
        rest = size2 % size1;
        StitchPoints(size2, step, rest, offSet1, offSet2, fitGeneralOffset, fitTotalSize, rotationFit, true);
    }
}

cv::Point3d GenerateSurfaceOld::VtkToCV(const double P1[3]) const
{
    cv::Point3d point = cv::Point3d(P1[0], P1[1], P1[2]);
    return point;
}

int GenerateSurfaceOld::GetDivisionAmount(vtkIdType size1, vtkIdType size2) const
{
    vtkIdType smollSize;
    if (size1 < size2)
    {
        smollSize = size1;
    }
    else
    {
        smollSize = size2;
    }
    if (smollSize >= 200)
    {
        return smollSize / 20;
    }
    else if (smollSize >= 60)
    {
        return smollSize / 10;
    }
    else if (smollSize >= 30)
    {
        return smollSize / 5;
    }
    else if (smollSize >= 21)
    {
        return 3;
    }
    else
    {
        return 2;
    }

}

vtkIdType GenerateSurfaceOld::FindByDistance(const cv::Point3d& currentFitPoint, const cv::Point3d& beforeFitPoint, const vtkSmartPointer<vtkPoints> Points, vtkIdType lastPos, vtkIdType step, double angleThreshold) const
{
    vtkIdType size1 = Points->GetNumberOfPoints();
    vtkIdType begin, end, saveEnd;

    cv::Point3d fitVector = beforeFitPoint - currentFitPoint;
    cv::Point3d fitVectorR = currentFitPoint - beforeFitPoint;
    cv::Point3d lastPoint = VtkToCV(Points->GetPoint(lastPos));
    cv::Point3d lastPoint2 = VtkToCV(Points->GetPoint(lastPos + 1));

    cv::Point3d initVector = lastPoint - currentFitPoint;

    begin = lastPos + 1;
    end = lastPos + 20 * step;
    saveEnd = lastPos + 10 * step;

    if (end >= size1)
    {
        end = size1 - 1;
    }

    if (saveEnd >= size1)
    {
        saveEnd = size1 - 1;
    }

    cv::Point3d changePoint, diff, changeVector;
    double distance;
    double distanceTemp = 9999999.0;
    double distanceTemp2 = 9999999.0;
    vtkIdType result = -1, result2;
    cv::Point3d rotVector;
    double angle, angle2;
    bool goodArea = false;
    for (vtkIdType i = begin; i <= end; i++)
    {
        double pnt[3];
        Points->GetPoint(i, pnt);
        changePoint = cv::Point3d(pnt[0], pnt[1], pnt[2]);

        rotVector = changePoint - currentFitPoint;

        diff = changePoint - currentFitPoint;
        distance = (diff.dot(diff));

        angle = GetAngleBetweenVectors(rotVector, fitVector) * 180.0 / vtkMath::Pi();

        angle2 = angle;

        angle = (abs(angle - 90.0));

        if (angleThreshold < 0)
        {
            if (angle < 20)
            {
                goodArea = true;
                if (distance < distanceTemp)
                {
                    distanceTemp = distance;
                    result = i;
                }
            }
            else if (goodArea == true)
            {
                break;
            }
        }
        else
        {
            if (angle2 >= angleThreshold )
            {
                result = i;
                break;
            }
        }

        if (distance < distanceTemp2 && i <= saveEnd)
        {
            distanceTemp2 = distance;
            result2 = i;
        }

    }

    if (result == -1)
    {
        result = result2;
    }
    return result;
}

void GenerateSurfaceOld::FindBestStepCombination(const vtkSmartPointer<vtkPoints> fitPoints, const vtkSmartPointer<vtkPoints> rotatePoints, std::vector<int>& fitSteps, std::vector<int>& rotateSteps) const
{
    if (fitPoints->GetNumberOfPoints() > rotatePoints->GetNumberOfPoints())
    {
        FindSteps(fitPoints, rotatePoints, fitSteps, rotateSteps);
    }
    else
    {
        FindSteps(rotatePoints, fitPoints, rotateSteps, fitSteps);
    }

}

void GenerateSurfaceOld::FindSteps(const vtkSmartPointer<vtkPoints> bigPoints, const vtkSmartPointer<vtkPoints> smollPoints, std::vector<int>& bigSteps, std::vector<int>& smollSteps) const
{
    vtkIdType size1 = smollPoints->GetNumberOfPoints();
    vtkIdType size2 = bigPoints->GetNumberOfPoints();

    int div = GetDivisionAmount(size1, size2);

    int diff = 9999999;

    int step1 = size1 / div;
    int step2 = size2 / div;

    int stepSize1, stepSize2;
    stepSize1 = step1;
    stepSize2 = 0;

    cv::Point3d fitPoint, fitA, fitB, beforeFitPoint, afterFitPoint;

    double angle = -1.0;
    while (div > 0)
    {
        angle = -1.0;
        if (div > 1)
        {
            fitPoint = VtkToCV(smollPoints->GetPoint(stepSize1 - 1));

            beforeFitPoint = VtkToCV(smollPoints->GetPoint(stepSize1 - step1));

            if (stepSize1 + step1 - 1 < size1 && size1 > 100 && size2 > 100)
            {
                afterFitPoint = VtkToCV(smollPoints->GetPoint(stepSize1 + step1 - 1));

                cv::Point3d beforeVector, nextVector;
                beforeVector = beforeFitPoint - fitPoint;
                nextVector = afterFitPoint - fitPoint;
                if ((GetAngleBetweenVectors(beforeVector, nextVector)) < closeAngle)
                {

                    angle = (GetAngleBetweenVectors(beforeVector, nextVector)) * 180.0 / vtkMath::Pi();
                    angle = (360.0 - angle) / 2.0;
                }
            }

            stepSize2 = FindByDistance(fitPoint, beforeFitPoint, bigPoints, stepSize2, step2, angle) + 1;

            if (stepSize2 > size2 - 3 || stepSize1 > size1 - 2)
            {
                stepSize2 = size2;
                stepSize1 = size1;
            }

        }

        smollSteps.push_back(stepSize1);
        bigSteps.push_back(stepSize2);

        div--;

        if (stepSize2 == size2)
        {
            break;
        }

        if (div == 1 || div == 0)
        {
            stepSize1 = size1;
            stepSize2 = size2;

            smollSteps.push_back(stepSize1);
            bigSteps.push_back(stepSize2);

            break;
        }
        else
        {
            stepSize1 = stepSize1 + step1;
            if (stepSize1 > size1)
            {
                stepSize1 = size1;
            }
        }
    }
}

vtkSmartPointer<vtkPoints> GenerateSurfaceOld::PrepareStitchTwoSlices(const vtkSmartPointer<vtkPoints> fitPoints, const vtkSmartPointer<vtkPoints> rotatePoints, vtkSmartPointer<vtkPoints>& outFitChange, bool remove)
{
    vtkIdType size1 = fitPoints->GetNumberOfPoints();
    vtkIdType size2 = rotatePoints->GetNumberOfPoints();
    vtkIdType generalOffset = allPoints->GetNumberOfPoints();
    vtkIdType offSet1, offSet2;
    vtkNew<vtkPoints> tempRotate;

    if (size1 == 0 || size2 == 0)
    {
        throw SegmentationException("There are empty slices. There can be no empty slices.");
    }

    vtkIdType rotFit;
    vtkIdType nearId = FindClosestPoint(fitPoints, rotatePoints, rotFit);

    for (vtkIdType i = nearId; i < size2; i++)
    {
        tempRotate->InsertNextPoint(rotatePoints->GetPoint(i));
    }

    for (vtkIdType i = 0; i < nearId; i++)
    {
        tempRotate->InsertNextPoint(rotatePoints->GetPoint(i));
    }

    for (vtkIdType i = rotFit; i < size1; i++)
    {
        outFitChange->InsertNextPoint(fitPoints->GetPoint(i));
    }

    for (vtkIdType i = 0; i < rotFit; i++)
    {
        outFitChange->InsertNextPoint(fitPoints->GetPoint(i));
    }

    int step1, stepSize1, stepSize2, fit, rot;
    int initFit = 0;
    int initRot = 0;
    vtkIdType beginFit, endFit, beginRot, endRot;

    std::vector<int> fitSteps, rotateSteps;

    FindBestStepCombination(outFitChange, tempRotate, fitSteps, rotateSteps);

    offSet1 = generalOffset;
    offSet2 = generalOffset + size1;
    for (int k = 0; k < fitSteps.size(); k++)
    {
        stepSize1 = fitSteps[k];
        stepSize2 = rotateSteps[k];

        fit = stepSize1 - initFit;
        rot = stepSize2 - initRot;

        MakeStitchTwoSlices(fit, rot, offSet1, offSet2);

        initFit = stepSize1 - 1;
        initRot = stepSize2 - 1;
        offSet1 = generalOffset + stepSize1 - 1;
        offSet2 = generalOffset + size1 + stepSize2 - 1;

    }

    for (vtkIdType i = 0; i < size1; i++)
    {
        allPoints->InsertNextPoint(outFitChange->GetPoint(i));
    }

    for (vtkIdType i = 0; i < size2; i++)
    {
        allPoints->InsertNextPoint(tempRotate->GetPoint(i));
    }

    beginFit = generalOffset;
    endFit = generalOffset + size1 - 1;

    beginRot = generalOffset + size1;
    endRot = generalOffset + size1 + size2 - 1;

    CreateCells(beginFit, endFit, beginRot);
    CreateCells(endFit, endRot, beginRot);
    return tempRotate;
}

void GenerateSurfaceOld::MakeStitchTwoSlices(vtkIdType size1, vtkIdType size2, vtkIdType offSet1, vtkIdType offSet2)
{
    int step, rest;

    if (size1 > size2)
    {
        step = size1 / size2;
        rest = size1 % size2;
        StitchPoints(size1, step, rest, offSet1, offSet2, 0, 0, 0);

    }
    else
    {
        step = size2 / size1;
        rest = size2 % size1;
        StitchPoints(size2, step, rest, offSet1, offSet2, 0, 0, 0, true);
    }
}

void GenerateSurfaceOld::MakeSurfaceFromContour(const vtkSmartPointer<vtkPoints> points, vtkIdType offset)
{
    vtkIdType size = points->GetNumberOfPoints();
    vtkNew<vtkCellArray> cells;
    for (unsigned int i = 0; i < size - 1; i++)
    {
        cells->InsertNextCell(2);
        cells->InsertCellPoint(i);
        cells->InsertCellPoint(i + 1);
    }

    vtkNew<vtkPolyData> contour;
    contour->SetPoints(points);
    contour->SetLines(cells);

    vtkNew<vtkContourTriangulator> triangulator;

    triangulator->SetInputData(contour);
    triangulator->Update();
    auto surface = triangulator->GetOutput();
    vtkNew<vtkIdList> myList;
    vtkSmartPointer<vtkCellArray> newCells = surface->GetPolys();

    while (newCells->GetNextCell(myList))
    {
        vtkIdType id1 = myList->GetId(0) + offset;
        vtkIdType id2 = myList->GetId(1) + offset;
        vtkIdType id3 = myList->GetId(2) + offset;
        CreateCells(id1, id2, id3);
    }
}

vtkSmartPointer<vtkPolyData> GenerateSurfaceOld::TriangularContours(const vtkSmartPointer<vtkPoints> a, const vtkSmartPointer<vtkPoints> b)
{
    int offset;
    SPlane plane;
    vtkNew<vtkPolyData> poly;
    vtkNew<vtkCellArray> myCells1;
    vtkNew<vtkPoints> p1;

    plane = GetFitPlane(a);
    for (unsigned int k = 0; k < a->GetNumberOfPoints(); k++)
    {
        p1->InsertNextPoint(a->GetPoint(k));

        myCells1->InsertNextCell(2);
        if (k == a->GetNumberOfPoints() - 1)
        {
            myCells1->InsertCellPoint(k);
            myCells1->InsertCellPoint(0);
        }
        else
        {
            myCells1->InsertCellPoint(k);
            myCells1->InsertCellPoint(k + 1);
        }
    }

    offset = p1->GetNumberOfPoints();

    for (int k = 0; k < b->GetNumberOfPoints(); k++)
    {
        cv::Point3d proj = (plane.getProjectionPoint(b->GetPoint(k)));
        double pnt[3] = { proj.x, proj.y, proj.z };
        p1->InsertNextPoint(pnt);
    }

    for (unsigned int k = offset; k < p1->GetNumberOfPoints(); k++)
    {
        myCells1->InsertNextCell(2);
        if (k == p1->GetNumberOfPoints() - 1)
        {
            myCells1->InsertCellPoint(k);
            myCells1->InsertCellPoint(offset);
        }
        else
        {
            myCells1->InsertCellPoint(k);
            myCells1->InsertCellPoint(k + 1);
        }
    }

    poly->SetPoints(p1);
    poly->SetLines(myCells1);

    vtkNew<vtkContourTriangulator> triangulator;
    triangulator->SetInputData(poly);
    triangulator->Update();

    auto surface = triangulator->GetOutput();

    return surface;
}

vtkSmartPointer<vtkPolyData> GenerateSurfaceOld::MakeSurface()
{
    vtkNew<vtkPolyData> result;
    vtkSmartPointer<vtkPoints> rotate;
    vtkSmartPointer<vtkPoints> outFitChange = vtkSmartPointer<vtkPoints>::New();

    for (int i = 0; i < slices.size() - 1; i++)
    {
        vtkNew<vtkPoints> p1, p2;
        vtkSmartPointer<vtkPoints> rotateChange;

        if (i == 0)
        {
            ExtractPoints(slices[i], p1);
            ExtractPoints(slices[i + 1], p2);
            rotateChange = PrepareStitchTwoSlices(p1, p2, outFitChange);
            MakeSurfaceFromContour(outFitChange);
        }
        else
        {
            ExtractPoints(slices[i + 1], p2);
            rotateChange = PrepareStitchOneSlices(rotate, p2);
        }

        rotate = rotateChange;

    }

    MakeSurfaceFromContour(rotate, allPoints->GetNumberOfPoints() - rotate->GetNumberOfPoints());

    result->SetPoints(allPoints);
    result->SetPolys(cells);

    return result;
}

bool GenerateSurfaceOld::IsProjectable(const vtkSmartPointer<vtkPolyData> a, const vtkSmartPointer<vtkPolyData> b) const
{
    double bound_a[6], bound_b[6];
    a->GetBounds(bound_a);
    b->GetBounds(bound_b);

    bool a1 = false, a2 = false, a3 = false;
    bool b1 = false, b2 = false, b3 = false;

    double average = 0.99;

    int pos0, pos1;
    pos0 = 0;
    pos1 = 1;

    if (abs(bound_a[pos0] - bound_a[pos1]) > EPSILON && abs(bound_b[pos0] - bound_b[pos1]) > EPSILON)
    {
        if (abs(bound_a[pos0] - bound_a[pos1]) < abs(bound_b[pos0] - bound_b[pos1]))
        {
            if (bound_a[pos0] > bound_b[pos0] && bound_a[pos1] < bound_b[pos1] && abs(bound_a[pos0] - bound_b[pos0]) >= average && abs(bound_a[pos1] - bound_b[pos1]) >= average)
            {
                a1 = true;
            }

        }
        else
        {
            if (bound_a[pos0] < bound_b[pos0] && bound_a[pos1] > bound_b[pos1] && abs(bound_a[pos0] - bound_b[pos0]) >= average && abs(bound_a[pos1] - bound_b[pos1]) >= average)
            {
                b1 = true;
            }
        }
    }
    else
    {
        if (abs(bound_a[pos0] - bound_a[pos1]) <= EPSILON && abs(bound_b[pos0] - bound_b[pos1]) <= EPSILON)
        {
            a1 = true;
            b1 = true;
        }
    }

    pos0 = 2;
    pos1 = 3;

    if (abs(bound_a[pos0] - bound_a[pos1]) > EPSILON && abs(bound_b[pos0] - bound_b[pos1]) > EPSILON)
    {
        if (abs(bound_a[pos0] - bound_a[pos1]) < abs(bound_b[pos0] - bound_b[pos1]))
        {
            if (bound_a[pos0] > bound_b[pos0] && bound_a[pos1] < bound_b[pos1] && abs(bound_a[pos0] - bound_b[pos0]) >= average && abs(bound_a[pos1] - bound_b[pos1]) >= average)
            {
                a2 = true;
            }

        }
        else
        {
            if (bound_a[pos0] < bound_b[pos0] && bound_a[pos1] > bound_b[pos1] && abs(bound_a[pos0] - bound_b[pos0]) >= average && abs(bound_a[pos1] - bound_b[pos1]) >= average)
            {
                b2 = true;
            }
        }
    }
    else
    {
        if (abs(bound_a[pos0] - bound_a[pos1]) <= EPSILON && abs(bound_b[pos0] - bound_b[pos1]) <= EPSILON)
        {
            a2 = true;
            b2 = true;
        }
    }

    pos0 = 4;
    pos1 = 5;

    if (abs(bound_a[pos0] - bound_a[pos1]) > EPSILON && abs(bound_b[pos0] - bound_b[pos1]) > EPSILON)
    {
        if (abs(bound_a[pos0] - bound_a[pos1]) < abs(bound_b[pos0] - bound_b[pos1]))
        {
            if (bound_a[pos0] > bound_b[pos0] && bound_a[pos1] < bound_b[pos1] && abs(bound_a[pos0] - bound_b[pos0]) >= average && abs(bound_a[pos1] - bound_b[pos1]) >= average)
            {
                a3 = true;
            }

        }
        else
        {
            if (bound_a[pos0] < bound_b[pos0] && bound_a[pos1] > bound_b[pos1] && abs(bound_a[pos0] - bound_b[pos0]) >= average && abs(bound_a[pos1] - bound_b[pos1]) >= average)
            {
                b3 = true;
            }
        }
    }
    else
    {
        if (abs(bound_a[pos0] - bound_a[pos1]) <= EPSILON && abs(bound_b[pos0] - bound_b[pos1]) <= EPSILON)
        {
            a3 = true;
            b3 = true;
        }
    }

    if (a1 == true && a2 == true && a3 == true)
    {
        return true;
    }

    if (b1 == true && b2 == true && b3 == true)
    {
        return true;
    }

    return false;
}

double GenerateSurfaceOld::GetAngleBetweenVectors(const cv::Point3d& V1, const cv::Point3d& V2, const SPlane& plane) const
{
    cv::Point3d a, b;
    if (plane.GetIsInit() == true)
    {
        cv::Point3d crossProduct1 = V1.cross(plane.getNormalVector());
        cv::Point3d crossProduct2 = V2.cross(plane.getNormalVector());

        if (abs(crossProduct1.dot(crossProduct1)) < EPSILON || abs(crossProduct2.dot(crossProduct2)) < EPSILON)
        {
            a = V1;
            b = V2;
        }
        else
        {
            a = plane.getProjectionVector(V1);
            b = plane.getProjectionVector(V2);
        }
    }
    else
    {
        a = V1;
        b = V2;
    }

    double scalar = a.dot(b);
    double magnitude = sqrt((a.dot(a)) * (b.dot(b)));
    double tCos = scalar / magnitude;
    if (tCos <= -1.0)
    {
        return vtkMath::Pi();
    }
    else if (tCos >= 1.0)
    {
        return 0;
    }
    else
    {
        double result = acos(tCos);
        return result;
    }
}

bool GenerateSurfaceOld::AreCollinear(const double point1[3], const double point2[3], const double point3[3]) const
{
    cv::Point3d P1, P2, P3;
    P1 = VtkToCV(point1);
    P2 = VtkToCV(point2);
    P3 = VtkToCV(point3);

    cv::Point3d directVector1, directVector2;
    directVector1 = P2 - P1;
    directVector2 = P3 - P1;

    cv::Point3d crossProduct = directVector1.cross(directVector2);

    if (abs(crossProduct.dot(crossProduct)) < EPSILON)
    {
        return true;
    }
    return false;
}

