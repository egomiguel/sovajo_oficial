
#include <vtkNew.h>
#include <fstream>
#include "SliceBorder.hpp"
#include "SliceBorderPrivate.hpp"
#include <opencv2/calib3d/calib3d.hpp>
#include "GenerateSurface.hpp"
#include "SegmentationException.hpp"
#include "vtkCutter.h"
#include "vtkPlane.h"
#include "vtkFlyingEdges3D.h"
#include "vtkPolyDataConnectivityFilter.h"

#include "vtkNamedColors.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include "vtkKdTreePointLocator.h"
#include "vtkSmoothPolyDataFilter.h"

SliceBorder::SliceBorder()
{
    m_data = new SliceBorderPrivate();
}

SliceBorder::~SliceBorder()
{
    m_data = NULL;
    delete m_data;
}

void SliceBorder::ExtractSortLines(const vtkSmartPointer<vtkPolyData> polyData, std::list<std::pair<vtkIdType, vtkIdType>>& lines) const
{
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
}

std::vector<PointTypeITK> SliceBorder::GetMainPoints(const vtkSmartPointer<vtkPolyData>& slice, double curvatureMaxAngle) const
{
    std::vector<PointTypeITK> result;
   
    vtkSmartPointer<vtkPoints> points = slice->GetPoints();

    std::list<std::pair<vtkIdType, vtkIdType>> lines;

    std::vector<cv::Point3d> sortPoints;

    ExtractSortLines(slice, lines);

    auto itl1 = lines.begin();
    auto itl2 = lines.end();

    for (; itl1 != itl2; ++itl1)
    {  
        double p[3];
        points->GetPoint((*itl1).first, p);
        sortPoints.push_back(cv::Point3d(p[0], p[1], p[2]));
    }

    std::vector<cv::Point3d> splinePoints = m_data->InterpolateSpline(sortPoints, 60);

    std::vector<cv::Point3d> mainPoints = m_data->GetContourMainPoints(splinePoints, curvatureMaxAngle);

    auto it1 = mainPoints.begin();
    auto it2 = mainPoints.end();

    for (; it1 != it2; ++it1)
    {
        PointTypeITK pointITK;
        pointITK[0] = (*it1).x;
        pointITK[1] = (*it1).y;
        pointITK[2] = (*it1).z;
        result.push_back(pointITK);
    }
    /////////////////////////////////////////////////////////////////////////////////////

   /* std::ofstream outfilePointBorder("BorderSegmentation.txt");
    for (int i = 0; i < slice->GetNumberOfPoints(); i++)
    {
        const auto pnt = slice->GetPoint(i);
        outfilePointBorder << pnt[0] << " " << pnt[1] << " " << pnt[2] << "\n";
    }
    outfilePointBorder.close();

    std::ofstream outfilePointSpline("BorderSegmentationSpline.txt");
    for (int i = 0; i < splinePoints.size(); i++)
    {
        cv::Point3d pnt = splinePoints[i];
        outfilePointSpline << pnt.x << " " << pnt.y << " " << pnt.z << "\n";
    }
    outfilePointSpline.close();

    std::ofstream outfilePoint("BorderSegmentationMainPoints.txt");
    for (int i = 0; i < mainPoints.size(); i++)
    {
        cv::Point3d pnt = mainPoints[i];
        outfilePoint << pnt.x << " " << pnt.y << " " << pnt.z << "\n";
    }
    outfilePoint.close();
    std::cout << "****Finish!!" << std::endl;*/
    /////////////////////////////////////////////////////////////////////////////////////////
    return result;
}

void SliceBorder::AddImproveSlice(const vtkSmartPointer<vtkPolyData>& slice)
{
    slices.push_back(slice);
}

vtkSmartPointer<vtkPolyData> SliceBorder::GetSurface3D(bool& result, std::string& msg) const
{
    GenerateSurface object(slices);

    vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();

    result = object.MakeSurface(surface, msg);

    return surface;
}

vtkSmartPointer<vtkPolyData> SliceBorder::GetSurface3D(const std::vector<vtkSmartPointer<vtkPolyData>>& polySlices, bool& result, std::string& msg) const
{
    GenerateSurface object(polySlices);
    
    vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();

    result = object.MakeSurface(surface, msg);

    return surface;
}

vtkSmartPointer<vtkPolyData> SliceBorder::GetSurface3D(const std::vector<vtkSmartPointer<vtkPolyData>>& polySlices, const double normal[3], bool& result, std::string& msg) const
{
    GenerateSurface object(polySlices, normal);

    vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();

    result = object.MakeSurface(surface, msg);

    return surface;
}

vtkSmartPointer<vtkPolyData> SliceBorder::SmoothSurface(const vtkSmartPointer<vtkPolyData> surface)
{
    vtkNew<vtkSmoothPolyDataFilter> smoothFilter;
    smoothFilter->SetInputData(surface);
    smoothFilter->SetNumberOfIterations(10);
    smoothFilter->SetRelaxationFactor(0.1);
    smoothFilter->FeatureEdgeSmoothingOff();
    smoothFilter->BoundarySmoothingOn();
    smoothFilter->Update();
    return smoothFilter->GetOutput();
}

std::vector<vtkSmartPointer<vtkPolyData>> SliceBorder::MakeSlicesFromPolyData(const vtkSmartPointer<vtkPolyData>& polyData, int axis, double distance) const
{
    double bounds[6];
    polyData->GetBounds(bounds);
    double x, y, z, min, max;
    vtkNew<vtkPlane> plane;
    vtkNew<vtkCutter> cutter;
    cutter->SetInputData(polyData);
    std::vector<vtkSmartPointer<vtkPolyData>> allSurface;

    if (distance < 1)
    {
        distance = 1;
    }

    if (axis == 1)
    {
        y = bounds[2];
        z = bounds[4];
        min = bounds[0];
        max = bounds[1];

        plane->SetNormal(1, 0, 0);

        double i;

        for (i = min + 0.2; i < max; i += distance)
        {
            plane->SetOrigin(i, y, z);
            cutter->SetCutFunction(plane);
            cutter->Update();

            auto contour = cutter->GetOutput();
            
            if (contour->GetNumberOfPoints() > 0)
            {
                vtkNew<vtkPolyDataConnectivityFilter> connectivityFilter;
                connectivityFilter->SetInputData(contour);
                connectivityFilter->SetExtractionModeToLargestRegion();
                connectivityFilter->Update();
                allSurface.push_back(connectivityFilter->GetOutput());
            }

        }

        if (i >= max)
        {
            if (max - (i - distance) >= 0.3)
            {
                i = max - 0.2;
                plane->SetOrigin(i, y, z);
                cutter->SetCutFunction(plane);
                cutter->Update();
                auto contour = cutter->GetOutput();

                if (contour->GetNumberOfPoints() > 0)
                {
                    vtkNew<vtkPolyDataConnectivityFilter> connectivityFilter;
                    connectivityFilter->SetInputData(contour);
                    connectivityFilter->SetExtractionModeToLargestRegion();
                    connectivityFilter->Update();
                    allSurface.push_back(connectivityFilter->GetOutput());
                }
            }
        }

        /*high = plane->EvaluateFunction(max, y, z);
        low = plane->EvaluateFunction(min, y, z);
        cutNumber = abs(max - min) / distance + 1;*/
    }
    else if (axis == 2)
    {
        x = bounds[0];
        z = bounds[4];
        min = bounds[2];
        max = bounds[3];

        plane->SetNormal(0, 1, 0);
        double i;

        for (i = min + 0.2; i < max; i += distance)
        {
            plane->SetOrigin(x, i, z);
            cutter->SetCutFunction(plane);
            cutter->Update();
            auto contour = cutter->GetOutput();

            if (contour->GetNumberOfPoints() > 0)
            {
                vtkNew<vtkPolyDataConnectivityFilter> connectivityFilter;
                connectivityFilter->SetInputData(contour);
                connectivityFilter->SetExtractionModeToLargestRegion();
                connectivityFilter->Update();
                allSurface.push_back(connectivityFilter->GetOutput());
            }
        }

        if (i >= max)
        {
            if (max - (i - distance) >= 0.3)
            {
                i = max - 0.2;

                plane->SetOrigin(x, i, z);
                cutter->SetCutFunction(plane);
                cutter->Update();
                auto contour = cutter->GetOutput();

                if (contour->GetNumberOfPoints() > 0)
                {
                    vtkNew<vtkPolyDataConnectivityFilter> connectivityFilter;
                    connectivityFilter->SetInputData(contour);
                    connectivityFilter->SetExtractionModeToLargestRegion();
                    connectivityFilter->Update();
                    allSurface.push_back(connectivityFilter->GetOutput());
                }
            }
        }

        /* high = plane->EvaluateFunction(x, max, z);
         low = plane->EvaluateFunction(x, min, z);
         cutNumber = abs(max - min) / distance + 1;*/
    }
    else
    {
        y = bounds[2];
        x = bounds[0];
        min = bounds[4];
        max = bounds[5];

        plane->SetNormal(0, 0, 1);
        double i;

        for (i = min + 0.2; i < max; i += distance)
        {
            plane->SetOrigin(x, y, i);
            cutter->SetCutFunction(plane);
            cutter->Update();
            auto contour = cutter->GetOutput();

            if (contour->GetNumberOfPoints() > 0)
            {
                vtkNew<vtkPolyDataConnectivityFilter> connectivityFilter;
                connectivityFilter->SetInputData(contour);
                connectivityFilter->SetExtractionModeToLargestRegion();
                connectivityFilter->Update();
                allSurface.push_back(connectivityFilter->GetOutput());
            }
        }

        if (i >= max)
        {
            if (max - (i - distance) >= 0.3)
            {
                i = max - 0.2;

                plane->SetOrigin(x, y, i);
                cutter->SetCutFunction(plane);
                cutter->Update();
                auto contour = cutter->GetOutput();

                if (contour->GetNumberOfPoints() > 0)
                {
                    vtkNew<vtkPolyDataConnectivityFilter> connectivityFilter;
                    connectivityFilter->SetInputData(contour);
                    connectivityFilter->SetExtractionModeToLargestRegion();
                    connectivityFilter->Update();
                    allSurface.push_back(connectivityFilter->GetOutput());
                }
            }
        }

        /*high = plane->EvaluateFunction(x, y, max);
        low = plane->EvaluateFunction(x, y, min);
        cutNumber = abs(max - min) / distance + 1;*/
    }

    //GenerateSurface object(allSurface);
    //auto result = object.MakeSurface();

    return allSurface;
}

//It is not finished
bool SliceBorder::MergeTwoPolyData(const vtkSmartPointer<vtkPolyData> poly1, const vtkSmartPointer<vtkPolyData> poly2, double normal[3], vtkSmartPointer<vtkPolyData>& result) const
{
    vtkNew<vtkKdTreePointLocator> kDTree;
    kDTree->SetDataSet(poly1);
    kDTree->BuildLocator();

    vtkSmartPointer<vtkPoints> pointsPoly2 = poly2->GetPoints();
    vtkSmartPointer<vtkPoints> pointsPoly1 = poly1->GetPoints();

    vtkIdType size2 = pointsPoly2->GetNumberOfPoints();

    vtkIdType id1, id2, tempId;
    double distance, tempDistance;
    tempDistance = 9999999.0;

    for (vtkIdType i = 0; i < size2; i++)
    {
        double pnt[3];
        pointsPoly2->GetPoint(i, pnt);
        tempId = kDTree->FindClosestPoint(pnt);

        distance = GetDistance(pointsPoly1->GetPoint(tempId), pnt);
        if (distance < tempDistance)
        {
            tempDistance = distance;
            id1 = tempId;
            id2 = i;
        }
    }

    poly1->BuildLinks();
    poly2->BuildLinks();

    vtkNew<vtkIdList> myList1, myList2;

    poly1->GetPointCells(id1, myList1);
    poly2->GetPointCells(id2, myList2);

    if (myList1->GetNumberOfIds() != 3 || myList2->GetNumberOfIds() != 3)
    {
        return false;
    }

    double pnt1[3];
    double pnt2[3];

    poly1->GetPoint(id1, pnt1);
    poly2->GetPoint(id2, pnt2);

    cv::Point3d normalPlane = m_data->ArrayToPoint(normal);
    cv::Point3d pointPlane1 = m_data->ArrayToPoint(pnt1);
    cv::Point3d pointPlane2 = m_data->ArrayToPoint(pnt1);

    return false;

}

//It is not finished
vtkSmartPointer<vtkPolyData> SliceBorder::RemoveBranchOnSlice(vtkSmartPointer<vtkPolyData> poly, double normal[3])
{
    vtkNew<vtkPolyData> result;

    vtkNew<vtkPolyDataConnectivityFilter> connectivityFilter;
    connectivityFilter->SetInputData(poly);
    connectivityFilter->SetExtractionModeToAllRegions();
    connectivityFilter->Update();

    int regCount = connectivityFilter->GetNumberOfExtractedRegions();

    if (regCount == 1)
    {
        return poly;
    }

    std::vector<vtkSmartPointer<vtkPolyData>> regions, regionsSort;

    for (int j = 0; j < regCount; j++)
    {
        connectivityFilter->SetInputData(poly);
        connectivityFilter->SetExtractionModeToSpecifiedRegions();
        connectivityFilter->InitializeSpecifiedRegionList();
        connectivityFilter->AddSpecifiedRegion(j);
        connectivityFilter->ScalarConnectivityOff();
        connectivityFilter->Update();
        regions.push_back(connectivityFilter->GetOutput());
    }

    double distance, tDistance;
    regionsSort.push_back(regions[0]);

    for (int i = 0; i < regions.size() - 1; i++)
    {
        tDistance = 99999999.0;
        double mainCenter[3];
        regions[i]->GetCenter(mainCenter);

        for (int j = i + 1; j < regions.size(); j++)
        {
            double center[3];
            regions[j]->GetCenter(center);

            distance = GetDistance(mainCenter, center);
            if (distance < tDistance)
            {
                tDistance = distance;
                regionsSort[i + 1] = regions[j];
            }
        }
    }
    return result;
}

double SliceBorder::GetDistance(const double p1[3], const double p2[3]) const
{
    double a = abs(p1[0] - p2[0]);
    double b = abs(p1[1] - p2[1]);
    double c = abs(p1[2] - p2[2]);

    return a * a + b * b + c * c;
}

std::vector<vtkSmartPointer<vtkPolyData>> SliceBorder::MakeSlicesFromPolyData(const vtkSmartPointer<vtkPolyData>& polyData, const SPlane& plane1, const cv::Point3d& pPoint, double distance) const
{
    cv::Point3d min = plane1.getProjectionPoint(pPoint);
    cv::Point3d max = pPoint;
    cv::Point3d vector = max - min;

    vector = vector / sqrt(vector.dot(vector));

    cv::Point3d diff = max - min;
    double lenght = sqrt(diff.dot(diff));

    vtkNew<vtkPlane> plane;
    vtkNew<vtkCutter> cutter;
    cutter->SetInputData(polyData);
    std::vector<vtkSmartPointer<vtkPolyData>> allSurface;

    if (distance < 1)
    {
        distance = 1;
    }

    plane->SetNormal(vector.x, vector.y, vector.z);

    double i = 0.2;
    cv::Point3d temp;

    while (i < lenght)
    {
        temp = min + i * vector;
        plane->SetOrigin(temp.x, temp.y, temp.z);
        cutter->SetCutFunction(plane);
        cutter->Update();

        auto contour = cutter->GetOutput();
        if (contour->GetNumberOfPoints() > 0)
        {
            vtkNew<vtkPolyDataConnectivityFilter> connectivityFilter;
            connectivityFilter->SetInputData(contour);
            connectivityFilter->SetExtractionModeToLargestRegion();
            connectivityFilter->Update();
            allSurface.push_back(connectivityFilter->GetOutput());
        }

        i = i + distance;
    }

    if (i >= lenght)
    {
        if (lenght - (i - distance) >= 0.3)
        {
            i = lenght - 0.2;

            temp = min + i * vector;
            plane->SetOrigin(temp.x, temp.y, temp.z);
            cutter->SetCutFunction(plane);
            cutter->Update();

            auto contour = cutter->GetOutput();
            if (contour->GetNumberOfPoints() > 0)
            {
                vtkNew<vtkPolyDataConnectivityFilter> connectivityFilter;
                connectivityFilter->SetInputData(contour);
                connectivityFilter->SetExtractionModeToLargestRegion();
                connectivityFilter->Update();
                allSurface.push_back(connectivityFilter->GetOutput());
            }
        }
    }

    return allSurface;

}

std::vector<vtkSmartPointer<vtkPolyData>> SliceBorder::MakeSlicesFromPolyData(const vtkSmartPointer<vtkPolyData>& polyData, const SPlane& plane1, const SPlane& plane2, double distance) const
{
    if (plane1.isParallelToPlane(plane2) == false)
    {
        throw SegmentationExceptionCode::PLANES_FOR_MAKE_SLICES_ARE_NOT_PARALLEL;
    }

    cv::Point3d min = plane1.getPoint();
    cv::Point3d max = plane2.getProjectionPoint(min);
    cv::Point3d vector = max - min;

    vector = vector / sqrt(vector.dot(vector));

    cv::Point3d diff = max - min;
    double lenght = sqrt(diff.dot(diff));

    vtkNew<vtkPlane> plane;
    vtkNew<vtkCutter> cutter;
    cutter->SetInputData(polyData);
    std::vector<vtkSmartPointer<vtkPolyData>> allSurface;

    if (distance < 1)
    {
        distance = 1;
    }

    plane->SetNormal(vector.x, vector.y, vector.z);

    double i = 0.2;
    cv::Point3d temp;

    while (i < lenght)
    {
        temp = min + i * vector;
        plane->SetOrigin(temp.x, temp.y, temp.z);
        cutter->SetCutFunction(plane);
        cutter->Update();

        auto contour = cutter->GetOutput();
        if (contour->GetNumberOfPoints() > 0)
        {
            vtkNew<vtkPolyDataConnectivityFilter> connectivityFilter;
            connectivityFilter->SetInputData(contour);
            connectivityFilter->SetExtractionModeToLargestRegion();
            connectivityFilter->Update();
            allSurface.push_back(connectivityFilter->GetOutput());
        }

        i = i + distance;
    }

    if (i >= lenght)
    {
        if (lenght - (i - distance) >= 0.3)
        {
            i = lenght - 0.2;

            temp = min + i * vector;
            plane->SetOrigin(temp.x, temp.y, temp.z);
            cutter->SetCutFunction(plane);
            cutter->Update();

            auto contour = cutter->GetOutput();
            if (contour->GetNumberOfPoints() > 0)
            {
                vtkNew<vtkPolyDataConnectivityFilter> connectivityFilter;
                connectivityFilter->SetInputData(contour);
                connectivityFilter->SetExtractionModeToLargestRegion();
                connectivityFilter->Update();
                allSurface.push_back(connectivityFilter->GetOutput());
            }
        }
    }

    return allSurface;

}

vtkSmartPointer<vtkPolyData> SliceBorder::GetCenterSlice(vtkSmartPointer<vtkPolyData> polyData, int axis)
{
    vtkNew<vtkPlane> plane;
    vtkNew<vtkCutter> cutter;

    if (axis == 1)
    {
        plane->SetNormal(1, 0, 0);
    }
    else if (axis == 2)
    {
        plane->SetNormal(0, 1, 0);
    }
    else
    {
        plane->SetNormal(0, 0, 1);
    }

    plane->SetOrigin(polyData->GetCenter());
    cutter->SetCutFunction(plane);
    cutter->SetInputData(polyData);
    cutter->Update();

    return cutter->GetOutput();
}

vtkSmartPointer<vtkPolyData> SliceBorder::ReadPolyData(std::string name)
{
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(name.c_str());
    reader->Update();
    return reader->GetOutput();
}

void SliceBorder::ShowPolyData(vtkSmartPointer<vtkPolyData> poly1, vtkSmartPointer<vtkPolyData> poly2)
{
    vtkNew<vtkNamedColors> colors;

    vtkNew<vtkPolyDataMapper> contoursMapper;
    contoursMapper->SetInputData(poly1);
    contoursMapper->ScalarVisibilityOff();

    vtkNew<vtkActor> contoursActor;
    contoursActor->SetMapper(contoursMapper);
    contoursActor->GetProperty()->SetRepresentationToWireframe();
    contoursActor->GetProperty()->ShadingOff();
    contoursActor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

    vtkNew<vtkPolyDataMapper> surfaceMapper;
    surfaceMapper->SetInputData(poly2);
    surfaceMapper->ScalarVisibilityOff();

    vtkNew<vtkActor> surfaceActor;
    surfaceActor->SetMapper(surfaceMapper);
    surfaceActor->GetProperty()->SetRepresentationToWireframe();
    surfaceActor->GetProperty()->ShadingOff();
    surfaceActor->GetProperty()->SetColor(
        colors->GetColor3d("MistyRose").GetData());

    // Create two renderers side by side to show the contours and the surface
    // separately
    //
    std::cout << "Press 't' for trackball interaction" << std::endl;
    std::cout << "Press 'r' to reset the camera" << std::endl;
    std::cout << "Press 'w' for wireframe representation" << std::endl;
    std::cout << "Press 's' for surface representation" << std::endl;

    vtkNew<vtkRenderer> renderer1;
    renderer1->SetViewport(0., 0., 0.5, 1.);
    renderer1->SetBackground(colors->GetColor3d("CadetBlue").GetData());

    vtkNew<vtkRenderer> renderer2;
    renderer2->SetViewport(0.5, 0., 1., 1.);
    renderer2->SetBackground(colors->GetColor3d("BurlyWood").GetData());

    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->SetSize(800, 400);
    renderWindow->SetWindowName("ContoursToSurface");

    renderWindow->AddRenderer(renderer1);
    renderWindow->AddRenderer(renderer2);

    vtkNew<vtkRenderWindowInteractor> interactor;
    interactor->SetRenderWindow(renderWindow);

    renderer1->AddViewProp(surfaceActor);
    renderer2->AddViewProp(contoursActor);
    renderWindow->Render();

    interactor->Start();
}

void SliceBorder::SavePolyData(vtkSmartPointer<vtkPolyData> poly1, std::string name)
{
    //vtk extention
    //vtkPolyData *myvtkPolyData = vtkPolyData::SafeDownCast(surfaceMapper->GetInputAsDataSet());
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(name.c_str());
    writer->SetInputData(poly1);
    writer->Write();
}