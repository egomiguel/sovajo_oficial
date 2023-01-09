#include "GenerateSurface.hpp"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkMath.h"
#include "Tiling.hpp"
#include "vtkContourTriangulator.h"
#include "vtkTriangle.h"


GenerateSurface::GenerateSurface(const std::vector<vtkSmartPointer<vtkPolyData>>& pSlices, const double* normal)
{
    if (normal)
    {
        if (normal[0] == 0 && normal[1] == 0 && normal[2] == 0)
        {
            EstimateNormal(pSlices);
        }
        else
        {
            double norm = sqrt((normal[0] * normal[0]) + (normal[1] * normal[1]) + (normal[2] * normal[2]));

            slicesNormal[0] = normal[0] / norm;
            slicesNormal[1] = normal[1] / norm;
            slicesNormal[2] = normal[2] / norm;
        }
    }
    else
    {
        EstimateNormal(pSlices);
    }

    auto it1 = pSlices.begin();
    auto it2 = pSlices.end();
    double vectorZ[3] = { 0, 0, 1 };
    for (; it1 != it2; ++it1)
    {
        zSlices.push_back(RotatePoly(*it1, slicesNormal, vectorZ));
    }
}

bool GenerateSurface::MakeSurface(vtkSmartPointer<vtkPolyData>& surface, std::string& msg)
{
    if (slicesNormal[0] == 0 && slicesNormal[1] == 0 && slicesNormal[2] == 0)
    {
        msg = "The normal vector to the slices cannot be estimated.";
        return false;
    }

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

    int offSet = 0;
    int size = zSlices.size();
    vtkSmartPointer<vtkPoints> tempUp;
    for (int i = 0; i < size - 1; i++)
    {
        if (i == 0)
        {
            Tiling myTiling(zSlices[i], zSlices[i + 1], offSet, true);
            myTiling.MakeTiling(points, cells);
            myTiling.GetTriangulateCellsDown(cells, offSet);
            tempUp = myTiling.GetOriginalUp();
        }
        else
        {
            Tiling myTiling(tempUp, zSlices[i + 1], offSet, true);
            myTiling.MakeTiling(points, cells, false);
            tempUp = myTiling.GetOriginalUp();

            if (i + 1 == size - 1)
            {
                offSet = points->GetNumberOfPoints() - tempUp->GetNumberOfPoints();
                myTiling.GetTriangulateCellsUp(cells, offSet);
            }
        }
        
        offSet = points->GetNumberOfPoints() - tempUp->GetNumberOfPoints();

    }
    vtkNew<vtkPolyData> resultRot;
    resultRot->SetPoints(points);
    resultRot->SetPolys(cells);

    double vectorZ[3] = { 0, 0, 1 };

    surface->DeepCopy(RotatePoly(resultRot, vectorZ, slicesNormal));

    return true;
}

double GenerateSurface::GetAngleBetweenVectors(const double v1[3], const double v2[3])
{
    double norm1 = sqrt((v1[0] * v1[0]) + (v1[1] * v1[1]) + (v1[2] * v1[2]));
    double norm2 = sqrt((v2[0] * v2[0]) + (v2[1] * v2[1]) + (v2[2] * v2[2]));

    double vector1[3] = { v1[0] / norm1, v1[1] / norm1, v1[2] / norm1 };
    double vector2[3] = { v2[0] / norm2, v2[1] / norm2, v2[2] / norm2 };

    double tCos = vtkMath::Dot(vector1, vector2);
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
        return acos(tCos);
    }
}

vtkSmartPointer<vtkPolyData> GenerateSurface::RotatePoly(vtkSmartPointer<vtkPolyData> poly, const double normal[3], const double newNormal[3])
{
    double theta = GetAngleBetweenVectors(normal, newNormal) * 180.0 / vtkMath::Pi();
    double axis[3] = { 0, 0, 0 };
    vtkMath::Cross(normal, newNormal, axis);

    vtkNew<vtkTransform> vtkTransform;
    vtkTransform->RotateWXYZ(theta, axis);

    vtkNew<vtkTransformPolyDataFilter> transformFilter;
    transformFilter->SetInputData(poly);
    transformFilter->SetTransform(vtkTransform);
    transformFilter->Update();

    auto resultTransform = transformFilter->GetOutput();
    return resultTransform;
}

void GenerateSurface::EstimateNormal(const std::vector<vtkSmartPointer<vtkPolyData>>& pSlices)
{
    if (pSlices.size() == 0)
    {
        slicesNormal[0] = 0;
        slicesNormal[1] = 0;
        slicesNormal[2] = 0;
        return;
    }

    auto it1 = pSlices.begin();
    auto it2 = pSlices.end();

    double x, y, z;

    std::vector<cv::Point3d> normalList;
    for (; it1 != it2; ++it1)
    {
        int size = (*it1)->GetNumberOfPoints();
        if (size > 2)
        {
            /*
            (*it1)->GetPoint(0, p1);
            (*it1)->GetPoint(size/2, p2);
            (*it1)->GetPoint(size - 1, p3);
            normalList.push_back(CalculateNormal(p1, p2, p3));*/

            double bound[6];
            (*it1)->GetBounds(bound);

            x = abs(bound[1] - bound[0]);
            y = abs(bound[3] - bound[2]);
            z = abs(bound[5] - bound[4]);

            double p1[3] = { bound[0], bound[2], bound[4] };
            double p2[3] = { bound[1], bound[3], bound[5] };

            if (x > y || x > z)
            {
                x = bound[0];
                y = bound[3];
                z = bound[5];
            }
            else if (y > x || y > z)
            {
                x = bound[1];
                y = bound[2];
                z = bound[5];
            }
            else
            {
                x = bound[1];
                y = bound[3];
                z = bound[4];
            }

            double p3[3] = { x, y, z };

            normalList.push_back(CalculateNormal(p1, p2, p3));
        }
        else
        {
            slicesNormal[0] = 0;
            slicesNormal[1] = 0;
            slicesNormal[2] = 0;
            return;
        }
    }

    cv::Point3d normalAmout = {0, 0, 0};
    double sum = 0;
    for (int i = 0; i < normalList.size() - 1; i++)
    {       
        cv::Point3d cross = (normalList[i]).cross(normalList[i + 1]);
        double amount = cross.dot(cross);

        if (amount < 0.0001 && amount > -0.0001)
        {
            normalAmout = normalAmout + normalList[i];
            sum++;
        }
        else
        {
            slicesNormal[0] = 0;
            slicesNormal[1] = 0;
            slicesNormal[2] = 0;
            return;
        }
    }

    normalAmout = normalAmout / sum;

    if (normalAmout.dot(normalAmout) == 0)
    {
        normalAmout = normalList[0];
    }

    normalAmout = normalAmout / sqrt(normalAmout.dot(normalAmout));

    slicesNormal[0] = normalAmout.x;
    slicesNormal[1] = normalAmout.y;
    slicesNormal[2] = normalAmout.z;

}

cv::Point3d GenerateSurface::ArrayToCv(const double pnt[3])
{
    cv::Point3d point = { pnt[0], pnt[1], pnt[2] };
    return point;
}

cv::Point3d GenerateSurface::CalculateNormal(const double pnt1[3], const double pnt2[3], const double pnt3[3])
{
    cv::Point3d p1 = ArrayToCv(pnt1);
    cv::Point3d p2 = ArrayToCv(pnt2);
    cv::Point3d p3 = ArrayToCv(pnt3);

    cv::Point3d v1 = p2 - p1;
    cv::Point3d v2 = p3 - p1;
    cv::Point3d v3 = v1.cross(v2);
    
    return v3/sqrt(v3.dot(v3));
}
