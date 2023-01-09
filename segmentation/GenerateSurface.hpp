#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vector"
#include <opencv2/calib3d/calib3d.hpp>

class GenerateSurface
{
private:
    double slicesNormal[3];
    std::vector<vtkSmartPointer<vtkPolyData>> zSlices;
    vtkSmartPointer<vtkPolyData> RotatePoly(vtkSmartPointer<vtkPolyData> poly, const double normal[3], const double newNormal[3]);
    void EstimateNormal(const std::vector<vtkSmartPointer<vtkPolyData>>& pSlices);
    cv::Point3d ArrayToCv(const double pnt[3]);
    cv::Point3d CalculateNormal(const double pnt1[3], const double pnt2[3], const double pnt3[3]);
    double GetAngleBetweenVectors(const double v1[3], const double v2[3]);
public:
    GenerateSurface(const std::vector<vtkSmartPointer<vtkPolyData>>& pSlices, const double* normal = nullptr);
    bool MakeSurface(vtkSmartPointer<vtkPolyData>& surface, std::string& msg);
};