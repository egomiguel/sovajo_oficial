#ifndef GENERATE_SURFACE_OLD_H
#define GENERATE_SURFACE_OLD_H

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vector"
#include "list"
#include "SPlane.hpp"
#include <opencv2/calib3d/calib3d.hpp>

class GenerateSurfaceOld
{
public:
    GenerateSurfaceOld(std::vector<vtkSmartPointer<vtkPolyData>> pSlices);
    vtkSmartPointer<vtkPolyData> MakeSurface();

private:
    std::vector<vtkSmartPointer<vtkPolyData>> slices;
    vtkSmartPointer<vtkPoints> allPoints;
    vtkSmartPointer<vtkCellArray> cells;

    SPlane GetFitPlane(const vtkSmartPointer<vtkPoints> pointsList) const;

    bool IsProjectable(const vtkSmartPointer<vtkPolyData> a, const vtkSmartPointer<vtkPolyData> b) const;

    vtkSmartPointer<vtkPolyData> TriangularContours(const vtkSmartPointer<vtkPoints> a, const vtkSmartPointer<vtkPoints> b);

    void FindBestStepCombination(const vtkSmartPointer<vtkPoints> fitPoints, const vtkSmartPointer<vtkPoints> rotatePoints, std::vector<int>& fitSteps, std::vector<int>& rotateSteps ) const;

    void FindSteps(const vtkSmartPointer<vtkPoints> bigPoints, const vtkSmartPointer<vtkPoints> smollPoints, std::vector<int>& bigSteps, std::vector<int>& smollSteps) const;


    cv::Point3d VtkToCV(const double P1[3]) const;
    vtkIdType FindByDistance(const cv::Point3d& currentFitPoint, const cv::Point3d& beforeFitPoint, const vtkSmartPointer<vtkPoints> Points, vtkIdType lastPos, vtkIdType step, double angleThreshold = -1.0) const;
    int GetDivisionAmount(vtkIdType size1, vtkIdType size2) const;
    bool AreCollinear(const double point1[3], const double point2[3], const double point3[3]) const;
    double GetAngleBetweenVectors(const cv::Point3d& V1, const cv::Point3d& V2, const SPlane& plane = SPlane()) const;
    double GetDistance(const double P1[3], const double P2[3], bool square = false) const;
    double GetDistance(const cv::Point3d& P1, const cv::Point3d& P2) const;
    double Resample(double lessDistance, const std::list<cv::Point3d>& linesIn, std::list<cv::Point3d>& linesOut) const;
    double ExtractSortLines(const vtkSmartPointer<vtkPolyData> polyData, std::list<std::pair<vtkIdType, vtkIdType>>& lines) const;
    double ExtractPoints(const vtkSmartPointer<vtkPolyData> polyData, vtkSmartPointer<vtkPoints> pPoints) const;
    std::list<cv::Point3d> ExtractPointsList(const vtkSmartPointer<vtkPolyData> polyData) const;
    void CreateCells(vtkIdType p1, vtkIdType p2, vtkIdType p3);
    vtkIdType FindClosestPoint(const vtkSmartPointer<vtkPoints> fitPoints, const vtkSmartPointer<vtkPoints> rotPoints, vtkIdType& rotFit) const;
    void StitchPoints(vtkIdType bigSetEnd, int step, int rest, vtkIdType generalOffSetBig, vtkIdType generalOffSetSmoll, vtkIdType fitGeneralOffset, vtkIdType fitSize, vtkIdType fitRot, bool localChangeOffset = false);
    vtkSmartPointer<vtkPoints> PrepareStitchTwoSlices(const vtkSmartPointer<vtkPoints> fitPoints, const vtkSmartPointer<vtkPoints> rotatePoints, vtkSmartPointer<vtkPoints>& outFitChange, bool remove = false);
    void MakeStitchTwoSlices(vtkIdType size1, vtkIdType size2, vtkIdType offSet1, vtkIdType offSet2);
    vtkSmartPointer<vtkPoints> PrepareStitchOneSlices(const vtkSmartPointer<vtkPoints> fitPoints, const vtkSmartPointer<vtkPoints> rotatePoints, bool remove = false);
    void MakeStitchOneSlices(vtkIdType size1, vtkIdType size2, vtkIdType offSet1, vtkIdType offSet2, vtkIdType fitGeneralOffset, vtkIdType fitTotalSize, vtkIdType rotationFit);
    void MakeSurfaceFromContour(const vtkSmartPointer<vtkPoints> points, vtkIdType offset = 0);
};

#endif