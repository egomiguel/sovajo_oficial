#include "Tiling.hpp"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkPointsProjectedHull.h"
#include "vtkKdTree.h"
#include "BinaryTree.hpp"
#include "vtkPolyLine.h"
#include "vtkLine.h"
#include "list"
#include "vtkKdTreePointLocator.h"
#include "vtkAppendPolyData.h"
#include "vtkContourTriangulator.h"


//void myShowPoint(const double a[3])
//{
//    std::cout << a[0] << "; " << a[1] << "; " << a[2] << std::endl;
//}
//
//void myShowPoint6(const double a[3])
//{
//    std::cout << a[0] << "; " << a[1] << "; " << a[2] << "; "<< a[3] << "; " << a[4] << "; " << a[5] << std::endl;
//}


Tiling::Tiling(const vtkSmartPointer<vtkPolyData> sliceDown, const vtkSmartPointer<vtkPolyData> sliceUp, int offSet, bool closed, bool makeNormalization)
{
    this->offSet = offSet;
    this->makeNormalization = makeNormalization;
    this->closed = closed;

    originalDown = ExtractSortPoints(sliceDown);
    originalUp = ExtractSortPoints(sliceUp);

    double downCenter[3];
    double upCenter[3];
    double downSize[6];
    double upSize[6];
    sliceDown->GetCenter(downCenter);
    sliceUp->GetCenter(upCenter);
    sliceDown->GetBounds(downSize);
    sliceUp->GetBounds(upSize);

    if (closed == true)
    {
        std::pair<int, int> rotationPair = GetNearestPoints(originalDown, originalUp);

        downRot = rotationPair.first;
        upRot = rotationPair.second;
    }
    else
    {
        downRot = 0;
        upRot = 0;
    }

    pointsDownRot = RotatePoints(originalDown, downRot);
    pointsUpRot = RotatePoints(originalUp, upRot);

    if (makeNormalization == true)
    {
        normaliceDownRot = normalice(pointsDownRot, downCenter, abs(downSize[1] - downSize[0]), abs(downSize[3] - downSize[2]));

        normaliceUpRot = normalice(pointsUpRot, upCenter, abs(upSize[1] - upSize[0]), abs(upSize[3] - upSize[2]));
    }
    else
    {
        normaliceDownRot = vtkSmartPointer<vtkPoints>::New();

        normaliceUpRot = vtkSmartPointer<vtkPoints>::New();
    }
}

Tiling::Tiling(const vtkSmartPointer<vtkPoints> pointsDown, const vtkSmartPointer<vtkPolyData> sliceUp, int offSet , bool closed, bool makeNormalization)
{
    this->offSet = offSet;
    this->makeNormalization = makeNormalization;
    this->closed = closed;

    originalDown = vtkSmartPointer<vtkPoints>::New();

    originalDown->DeepCopy(pointsDown);
    originalUp = ExtractSortPoints(sliceUp);

    vtkNew<vtkPolyData> zDown;
    zDown->SetPoints(originalDown);

    double downCenter[3];
    double upCenter[3];
    double downSize[6];
    double upSize[6];
    zDown->GetCenter(downCenter);
    sliceUp->GetCenter(upCenter);
    zDown->GetBounds(downSize);
    sliceUp->GetBounds(upSize);

    if (closed == true)
    {
        std::pair<int, int> rotationPair = GetNearestPoints(originalDown, originalUp);

        downRot = rotationPair.first;
        upRot = rotationPair.second;
    }
    else
    {
        downRot = 0;
        upRot = 0;
    }

    pointsDownRot = RotatePoints(originalDown, downRot);
    pointsUpRot = RotatePoints(originalUp, upRot);

    if (makeNormalization == true)
    {
        normaliceDownRot = normalice(pointsDownRot, downCenter, abs(downSize[1] - downSize[0]), abs(downSize[3] - downSize[2]));

        normaliceUpRot = normalice(pointsUpRot, upCenter, abs(upSize[1] - upSize[0]), abs(upSize[3] - upSize[2]));
    }
    else
    {
        normaliceDownRot = vtkSmartPointer<vtkPoints>::New();

        normaliceUpRot = vtkSmartPointer<vtkPoints>::New();
    }
}

Tiling::Tiling(const vtkSmartPointer<vtkPoints> pointsDown, const vtkSmartPointer<vtkPoints> pointsUp, int offSet, bool closed, bool makeNormalization)
{
    this->offSet = 0;
    this->makeNormalization = makeNormalization;
    this->closed = closed;

    originalDown = vtkSmartPointer<vtkPoints>::New();
    originalUp = vtkSmartPointer<vtkPoints>::New();

    originalDown->DeepCopy(pointsDown);
    originalUp->DeepCopy(pointsUp);

    vtkNew<vtkPolyData> zDown, zUp;
    zDown->SetPoints(originalDown);
    zUp->SetPoints(originalUp);

    double downCenter[3];
    double upCenter[3];
    double downSize[6];
    double upSize[6];
    zDown->GetCenter(downCenter);
    zUp->GetCenter(upCenter);
    zDown->GetBounds(downSize);
    zUp->GetBounds(upSize);

    if (closed == true)
    {
        std::pair<int, int> rotationPair = GetNearestPoints(originalDown, originalUp);

        downRot = rotationPair.first;
        upRot = rotationPair.second;
    }
    else
    {
        downRot = 0;
        upRot = 0;
    }

    pointsDownRot = RotatePoints(originalDown, downRot);
    pointsUpRot = RotatePoints(originalUp, upRot);

    normaliceDownRot = normalice(pointsDownRot, downCenter, abs(downSize[1] - downSize[0]), abs(downSize[3] - downSize[2]));

    normaliceUpRot = normalice(pointsUpRot, upCenter, abs(upSize[1] - upSize[0]), abs(upSize[3] - upSize[2]));
}

std::pair<int, int> Tiling::GetNearestPoints(const vtkSmartPointer<vtkPoints> pointsDown, const vtkSmartPointer<vtkPoints> pointsUp)
{
    vtkNew<vtkPolyData> polydata;
    polydata->SetPoints(pointsDown);

    vtkNew<vtkKdTreePointLocator> kDTree;
    kDTree->SetDataSet(polydata);
    kDTree->BuildLocator();

    vtkIdType size = pointsUp->GetNumberOfPoints();

    vtkIdType downId, upId, tempId;
    double distance, tempDistance;
    distance = 9999999.0;

    for (vtkIdType i = 0; i < size; i++)
    {
        double pnt[3];
        pointsUp->GetPoint(i, pnt);
        tempId = kDTree->FindClosestPoint(pnt);

        tempDistance = GetDistance(pointsDown->GetPoint(tempId), pnt);
        if (tempDistance < distance)
        {
            distance = tempDistance;
            downId = tempId;
            upId = i;
        }
    }
    
    return std::pair<int, int>(downId, upId);
}

vtkSmartPointer<vtkPoints> Tiling::normalice(vtkSmartPointer<vtkPoints> points, double center[3], double sizeX, double sizeY)
{
    vtkNew<vtkPoints> normalicePoints;
    double x, y;
    for (int i = 0; i < points->GetNumberOfPoints(); i++)
    {
        double pnt[3];
        points->GetPoint(i, pnt);

        x = (pnt[0] - center[0]) / sizeX;
        y = (pnt[1] - center[1]) / sizeY;
        normalicePoints->InsertNextPoint(x, y, pnt[2]);
    }
    return normalicePoints;
}

vtkSmartPointer<vtkPoints> Tiling::ExtractSortPoints(const vtkSmartPointer<vtkPolyData> polyData) const
{
    std::list<std::pair<vtkIdType, vtkIdType>> lines;
    vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
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

    vtkNew<vtkPoints> result;

    auto it1 = lines.begin();
    auto it2 = lines.end();

    for (; it1 != it2; ++it1)
    {
        double p1[3];
        points->GetPoint(it1->first, p1);
        result->InsertNextPoint(p1);
    }

    if (lines.size() > 0)
    {
        if (lines.front().first != lines.back().second)
        {
            result->InsertNextPoint(points->GetPoint(lines.back().second));
        }
    }

    return result;
}

std::vector<std::pair<int, int>> Tiling::GetTilingPoints(const vtkSmartPointer<vtkPoints> downPoints, const vtkSmartPointer<vtkPoints> upPoints)
{
    GenerateTree myTree(downPoints, upPoints, closed);
    std::vector<std::pair<int, int>> path = myTree.GetPath();

    auto it1 = path.begin();
    auto it2 = path.end();

    std::reverse(path.begin(), path.end());
    
    return path;
}

void Tiling::MakeTiling(vtkSmartPointer<vtkPoints>& points, vtkSmartPointer<vtkCellArray>& cells, bool addFixPoints)
{
    std::vector<std::pair<int, int>> path;

    if (makeNormalization == true)
    {
        path = GetTilingPoints(normaliceDownRot, normaliceUpRot);
    }
    else
    {
        path = GetTilingPoints(pointsDownRot, pointsUpRot);
    }

    UpdatePath(path);

    int offSetUp = originalDown->GetNumberOfPoints() + offSet;

    if (addFixPoints == true)
    {
        for (int i = 0; i < originalDown->GetNumberOfPoints(); i++)
        {
            points->InsertNextPoint(originalDown->GetPoint(i));
        }
    }

    for (int i = 0; i < originalUp->GetNumberOfPoints(); i++)
    {
        points->InsertNextPoint(originalUp->GetPoint(i));
    }

    auto it1 = path.begin();
    auto it2 = path.end();
    auto finalPos = std::next(it1, path.size() - 1);

    for (; it1 != finalPos; ++it1)
    {
        auto nextIt = std::next(it1, 1);
        cells->InsertNextCell(MakeCell(*it1, *nextIt, offSet, offSetUp));
    }

}


vtkSmartPointer<vtkPolyData> Tiling::MakeTiling()
{
    std::vector<std::pair<int, int>> path;

    if (makeNormalization == true)
    {
        path = GetTilingPoints(normaliceDownRot, normaliceUpRot);
    }
    else
    {
        path = GetTilingPoints(pointsDownRot, pointsUpRot);
    }

    UpdatePath(path);

    vtkNew<vtkPoints> allPoints;
    vtkNew<vtkCellArray> cells;

    int offSetUp = originalDown->GetNumberOfPoints() + offSet;

    for (int i = 0; i < originalDown->GetNumberOfPoints(); i++)
    {
        allPoints->InsertNextPoint(originalDown->GetPoint(i));
    }

    for (int i = 0; i < originalUp->GetNumberOfPoints(); i++)
    {
        allPoints->InsertNextPoint(originalUp->GetPoint(i));
    }

    vtkNew<vtkPolyData> result;
    result->SetPoints(allPoints);

    auto it1 = path.begin();
    auto it2 = path.end();
    auto finalPos = std::next(it1, path.size() - 1);

    for (; it1 != finalPos; ++it1)
    {
        auto nextIt = std::next(it1, 1);
        cells->InsertNextCell(MakeCell(*it1, *nextIt, offSet, offSetUp));
    }

    result->SetPolys(cells);
    return result;
}


vtkSmartPointer<vtkTriangle> Tiling::MakeCell(const std::pair<int, int>& node1, const std::pair<int, int>& node2, int offSetDown, int offSetUp)
{
    vtkNew<vtkTriangle> triangle;
    int a, b, c;

    if (node1.first != node2.first)
    {
        a = node1.first;
        b = node2.first;
        c = node1.second;

        a = a + offSetDown;
        b = b + offSetDown;
        c = c + offSetUp;
    }
    else
    {
        a = node1.second;
        b = node2.second;
        c = node1.first;

        a = a + offSetUp;
        b = b + offSetUp;
        c = c + offSetDown;
    }

    triangle->GetPointIds()->SetId(0, a);
    triangle->GetPointIds()->SetId(1, b);
    triangle->GetPointIds()->SetId(2, c);

    return triangle;
}

double Tiling::GetDistance(const double p1[3], const double p2[3]) const
{
    double a = abs(p1[0] - p2[0]);
    double b = abs(p1[1] - p2[1]);
    double c = abs(p1[2] - p2[2]);

    return a * a + b * b + c * c;
}

vtkSmartPointer<vtkPoints> Tiling::RotatePoints(vtkSmartPointer<vtkPoints> points, vtkIdType id)
{
    if (id == 0)
    {
        return points;
    }

    vtkNew<vtkPoints> result;
    vtkIdType size = points->GetNumberOfPoints();

    for (int i = id; i < size; i++)
    {
        result->InsertNextPoint(points->GetPoint(i));
    }

    for (int i = 0; i < id; i++)
    {
        result->InsertNextPoint(points->GetPoint(i));
    }

    return result;
}

void Tiling::UpdatePath(std::vector<std::pair<int, int>>& path)
{
    auto it1 = path.begin();
    auto it2 = path.end();

    int downSize = pointsDownRot->GetNumberOfPoints();
    int upSize = pointsUpRot->GetNumberOfPoints();

    for (; it1 != it2; ++it1)
    {
        int a = it1->first;
        int b = it1->second;

        a = a + downRot;
        b = b + upRot;

        if (a >= downSize)
        {
            a = a - downSize;
        }

        if (b >= upSize)
        {
            b = b - upSize;
        }

        it1->first = a;
        it1->second = b;
    }
}

vtkSmartPointer<vtkPolyData> Tiling::PointsToContour(const vtkSmartPointer<vtkPoints> points)
{
    vtkNew<vtkCellArray> cells;

    for (unsigned int i = 0; i < points->GetNumberOfPoints(); i++)
    {
        vtkNew<vtkLine> myLine;
        if (i == points->GetNumberOfPoints() - 1)
        {
            myLine->GetPointIds()->SetId(0, i);
            myLine->GetPointIds()->SetId(1, 0);
        }
        else
        {
            myLine->GetPointIds()->SetId(0, i);
            myLine->GetPointIds()->SetId(1, i + 1);
        }
        cells->InsertNextCell(myLine);
    }

    vtkNew<vtkPolyData> polyData;
    polyData->SetPoints(points);
    polyData->SetLines(cells);
    return polyData;
}

vtkSmartPointer<vtkPolyData> Tiling::TriangulateContour(const vtkSmartPointer<vtkPolyData> contour)
{
    vtkNew<vtkContourTriangulator> triangulator;
    triangulator->SetInputData(contour);
    triangulator->Update();
    auto surface = triangulator->GetOutput();
    return surface;
}

void Tiling::GetTriangulateCellsDown(vtkSmartPointer<vtkCellArray>& cells, int currentOffset)
{
    vtkSmartPointer<vtkPolyData> myConTour = PointsToContour(originalDown);
    vtkSmartPointer<vtkPolyData> contourTriangulation = TriangulateContour(myConTour);

    vtkNew<vtkIdList> myList;
    vtkSmartPointer<vtkCellArray> newCells = contourTriangulation->GetPolys();

    while (newCells->GetNextCell(myList))
    {
        vtkIdType id1 = myList->GetId(0) + currentOffset;
        vtkIdType id2 = myList->GetId(1) + currentOffset;
        vtkIdType id3 = myList->GetId(2) + currentOffset;
        cells->InsertNextCell(CreateCells(id1, id2, id3));
    }
}

void Tiling::GetTriangulateCellsUp(vtkSmartPointer<vtkCellArray>& cells, int currentOffset)
{
    vtkSmartPointer<vtkPolyData> myConTour = PointsToContour(originalUp);
    vtkSmartPointer<vtkPolyData> contourTriangulation = TriangulateContour(myConTour);

    vtkNew<vtkIdList> myList;
    vtkSmartPointer<vtkCellArray> newCells = contourTriangulation->GetPolys();

    while (newCells->GetNextCell(myList))
    {
        vtkIdType id1 = myList->GetId(0) + currentOffset;
        vtkIdType id2 = myList->GetId(1) + currentOffset;
        vtkIdType id3 = myList->GetId(2) + currentOffset;
        cells->InsertNextCell(CreateCells(id1, id2, id3));
    }
}

vtkSmartPointer<vtkTriangle> Tiling::CreateCells(vtkIdType p1, vtkIdType p2, vtkIdType p3)
{
    vtkNew<vtkTriangle> triangle;
    triangle->GetPointIds()->SetId(0, p1);
    triangle->GetPointIds()->SetId(1, p2);
    triangle->GetPointIds()->SetId(2, p3);
    return triangle;
}

vtkSmartPointer<vtkPoints> Tiling::GetOriginalUp() const
{
   /* vtkNew<vtkPoints> result;
    result->DeepCopy(originalUp);
    return result;*/
    return originalUp;
}

vtkSmartPointer<vtkPoints> Tiling::GetOriginalDown() const
{
    /*vtkNew<vtkPoints> result;
    result->DeepCopy(originalDown);
    return result;*/
    return originalDown;
}