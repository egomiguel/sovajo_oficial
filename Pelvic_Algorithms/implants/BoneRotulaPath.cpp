#include "BoneRotulaPath.hpp"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkSphereSource.h"
#include "vtkAppendPolyData.h"
#include "FitParabola.hpp"
#include "vtkImplicitPolyDataDistance.h"
#include <Eigen/Core>
#include "ImplantTools.hpp"
#include <unsupported/Eigen/Splines>
//#include <pcl/surface/on_nurbs/fitting_curve_pdm.h>
#include <fstream>

BoneRotulaPath::BoneRotulaPath(const Knee& knee)
{
    femurPoly = knee.GetFemurPoly();
    mKneeGroovePath = knee.getKneeGroovePath();
    //femurCoordenate = new CoordenateSystemFemur(knee.getHipCenter(), knee.getLateralEpicondyle(), knee.getMedialEpicondylePerp(), knee.getFemurKneeCenter());
}

BoneRotulaPath::~BoneRotulaPath()
{
    /*delete femurCoordenate;
    femurCoordenate = NULL;*/
}

std::vector<PointTypeITK> BoneRotulaPath::getKneeCapPathPoints() const
{
    std::vector<PointTypeITK> path;
    for (int i = 0; i < mKneeGroovePath.size(); i++)
    {
        path.push_back(mKneeGroovePath[i].ToITKPoint());
    }
    return path;
}

vtkSmartPointer<vtkPolyData> BoneRotulaPath::CreateSphereTest(const Point& pPoint)
{
    double pnt[3];
    pnt[0] = pPoint.x;
    pnt[1] = pPoint.y;
    pnt[2] = pPoint.z;

    vtkNew<vtkSphereSource> sphere;
    sphere->SetCenter(pnt);
    sphere->SetRadius(1);
    sphere->Update();

    return sphere->GetOutput();
}

vtkSmartPointer<vtkPolyData> BoneRotulaPath::JoinSagitalPointTest()
{

    vtkNew<vtkAppendPolyData> filter;
    filter->AddInputData(femurPoly);

    for (int i = 0; i < mKneeGroovePath.size(); i++)
    {
        filter->AddInputData(CreateSphereTest(mKneeGroovePath[i]));
    }

    filter->Update();

    return filter->GetOutput();
}


/*
void BoneRotulaPath::FitPoints()
{
    std::vector<Point> pPoints;
    pcl::on_nurbs::NurbsDataCurve curve_data;
    pcl::on_nurbs::vector_vec3d data;
    for (int i = 1; i <= 10; i++)
    {
        pPoints.push_back(Point(i, i*i, 0));
        pPoints.push_back(Point(i, (i-1)*(i-1), 0));

        Point P1(i, i*i, 0);
        Point P2(i, (i - 1)*(i - 1), 0);

        data.emplace_back(P1.x, P1.y, P1.z);
        data.emplace_back(P2.x, P2.y, P2.z);
    }

    RotulaSavePoints(pPoints, "rotula_fit_points.txt");

    curve_data.interior = data;

    unsigned order(3);
    unsigned n_control_points(20);

    pcl::on_nurbs::FittingCurve::Parameter curve_params;
    curve_params.smoothness = 0.000001;

    ON_NurbsCurve curve = pcl::on_nurbs::FittingCurve::initNurbsCurvePCA(order, curve_data.interior, n_control_points);

    pcl::on_nurbs::FittingCurve fit(&curve_data, curve);
    fit.assemble(curve_params);
    fit.solve();

    ON_NurbsCurve finalCurve = fit.m_nurbs;

    std::vector<Point> result;

    for (int i = 0; i < 20; i++)
    {
        ON_3dPoint cp;
        curve.GetCV(i, cp);
        result.push_back(Point(cp.x, cp.y, cp.z));
    }

    RotulaSavePoints(result, "rotula_fit_points.txt");*/

    /*RotulaSavePoints(pPoints, "rotula_init_points.txt");

    std::vector<cv::Point2d> result;

    int amount = pPoints.size();

    Eigen::MatrixXd points(2, pPoints.size());
    int row_index = 0;
    for (auto const way_point : pPoints)
    {
        points.col(row_index) << way_point.x, way_point.y;
        row_index++;
    }
    Eigen::Spline2d spline = Eigen::SplineFitting<Eigen::Spline2d>::Interpolate(points, 2);

    double time_ = 0;
    for (int i = 0; i < amount; i++) {
        time_ += 1.0 / (amount * 1.0);
        Eigen::VectorXd values = spline(time_);
        result.push_back(cv::Point2d(values[0], values[1]));
    }

    RotulaSavePoints(result, "rotula_fit_points.txt");
}
*/