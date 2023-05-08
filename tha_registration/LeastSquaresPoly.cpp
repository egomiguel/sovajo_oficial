#include "LeastSquaresPoly.hpp"

///////////////////////////////////////////////////////////////////////////////

#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>
#include <vtkNamedColors.h>
#include <vtkPolyDataMapper.h>

using namespace THA::RIGISTRATION;

/*
vtkSmartPointer<vtkPolyData> CreateSphereTest(const cv::Point3d& pPoint)
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

vtkSmartPointer<vtkPolyData> MergePolyWithSphere(const vtkSmartPointer<vtkPolyData> poly, const std::vector<cv::Point3d>& pPoints)
{
    vtkNew<vtkAppendPolyData> filter;
    filter->AddInputData(poly);

    for (int i = 0; i < pPoints.size(); i++)
    {
        vtkSmartPointer<vtkPolyData> mySphere = CreateSphereTest(pPoints[i]);
        filter->AddInputData(mySphere);
    }

    filter->Update();
    return filter->GetOutput();
}

void ShowPolyData(vtkSmartPointer<vtkPolyData> poly1, vtkSmartPointer<vtkPolyData> poly2)
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

void showProcess(vtkSmartPointer<vtkPolyData> poly, const std::vector<cv::Point3d>& points)
{
    vtkNew<vtkNamedColors> colors;

    vtkNew<vtkPolyDataMapper> contoursMapper;
    contoursMapper->SetInputData(poly);
    contoursMapper->ScalarVisibilityOff();

    vtkNew<vtkActor> polyActor;
    polyActor->SetMapper(contoursMapper);
    polyActor->GetProperty()->SetRepresentationToWireframe();
    polyActor->GetProperty()->ShadingOff();
    polyActor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

    std::vector<vtkSmartPointer<vtkActor>> pointsActor;

    for (int i = 0; i < points.size(); i++)
    {
        double pnt[3];
        pnt[0] = points[i].x;
        pnt[1] = points[i].y;
        pnt[2] = points[i].z;

        vtkNew<vtkSphereSource> sphere;
        sphere->SetCenter(pnt);
        sphere->SetRadius(1);
        sphere->Update();

        vtkNew<vtkPolyDataMapper> sphereMapper;
        sphereMapper->SetInputData(sphere->GetOutput());
        sphereMapper->ScalarVisibilityOff();

        vtkNew<vtkActor> sphereActor;
        sphereActor->SetMapper(sphereMapper);
        sphereActor->GetProperty()->SetRepresentationToWireframe();
        sphereActor->GetProperty()->ShadingOff();
        sphereActor->GetProperty()->SetColor(colors->GetColor3d("blue").GetData());

        pointsActor.push_back(sphereActor);
    }

    vtkNew<vtkRenderer> renderer;
    //renderer->SetViewport(0., 0., 0.5, 1.);
    renderer->SetBackground(colors->GetColor3d("CadetBlue").GetData());

    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->SetSize(800, 400);
    renderWindow->SetWindowName("Surface");

    renderWindow->AddRenderer(renderer);

    vtkNew<vtkRenderWindowInteractor> interactor;
    interactor->SetRenderWindow(renderWindow);

    renderer->AddActor(polyActor);

    for (int i = 0; i < pointsActor.size(); i++)
    {
        renderer->AddActor(pointsActor[i]);
    }

    renderWindow->Render();

    interactor->Start();
}
*/
////////////////////////////////////////////////////////////////////////////////////


LeastSquaresPoly::LeastSquaresPoly(const std::vector<PointTypeITK>& pSourcePoints, const vtkSmartPointer<vtkPolyData>& pSurface)
{
    for (int i = 0; i < pSourcePoints.size(); i++)
    {
        source.push_back(cv::Point3d(pSourcePoints[i][0], pSourcePoints[i][1], pSourcePoints[i][2]));
    }

    surface = pSurface;
    implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    implicitPolyDataDistance->SetInput(surface);
    maxError = 0.35;
}

LeastSquaresPoly::LeastSquaresPoly(const std::vector<cv::Point3d>& pSourcePoints, const vtkSmartPointer<vtkPolyData>& pSurface)
{
    source = pSourcePoints;
    surface = pSurface;
    implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    implicitPolyDataDistance->SetInput(surface);
    maxError = 0.35;
}

cv::Mat LeastSquaresPoly::Rx(double angle)
{
    cv::Mat matrix = cv::Mat::zeros(3, 3, CV_64F);

    matrix.at<double>(0, 0) = 1;
    matrix.at<double>(0, 1) = 0;
    matrix.at<double>(0, 2) = 0;

    matrix.at<double>(1, 0) = 0;
    matrix.at<double>(1, 1) = cos(angle);
    matrix.at<double>(1, 2) = -sin(angle);

    matrix.at<double>(2, 0) = 0;
    matrix.at<double>(2, 1) = sin(angle);
    matrix.at<double>(2, 2) = cos(angle);
    return matrix;
}

cv::Mat LeastSquaresPoly::Ry(double angle)
{
    cv::Mat matrix = cv::Mat::zeros(3, 3, CV_64F);

    matrix.at<double>(0, 0) = cos(angle);
    matrix.at<double>(0, 1) = 0;
    matrix.at<double>(0, 2) = sin(angle);

    matrix.at<double>(1, 0) = 0;
    matrix.at<double>(1, 1) = 1;
    matrix.at<double>(1, 2) = 0;

    matrix.at<double>(2, 0) = -sin(angle);
    matrix.at<double>(2, 1) = 0;
    matrix.at<double>(2, 2) = cos(angle);
    return matrix;
}

cv::Mat LeastSquaresPoly::Rz(double angle)
{
    cv::Mat matrix = cv::Mat::zeros(3, 3, CV_64F);

    matrix.at<double>(0, 0) = cos(angle);
    matrix.at<double>(0, 1) = -sin(angle);
    matrix.at<double>(0, 2) = 0;

    matrix.at<double>(1, 0) = sin(angle);
    matrix.at<double>(1, 1) = cos(angle);
    matrix.at<double>(1, 2) = 0;

    matrix.at<double>(2, 0) = 0;
    matrix.at<double>(2, 1) = 0;
    matrix.at<double>(2, 2) = 1;
    return matrix;
}

cv::Mat LeastSquaresPoly::DRx(double angle)
{
    cv::Mat matrix = cv::Mat::zeros(3, 3, CV_64F);

    matrix.at<double>(0, 0) = 0;
    matrix.at<double>(0, 1) = 0;
    matrix.at<double>(0, 2) = 0;

    matrix.at<double>(1, 0) = 0;
    matrix.at<double>(1, 1) = -sin(angle);
    matrix.at<double>(1, 2) = -cos(angle);

    matrix.at<double>(2, 0) = 0;
    matrix.at<double>(2, 1) = cos(angle);
    matrix.at<double>(2, 2) = -sin(angle);
    return matrix;
}

cv::Mat LeastSquaresPoly::DRy(double angle)
{
    cv::Mat matrix = cv::Mat::zeros(3, 3, CV_64F);

    matrix.at<double>(0, 0) = -sin(angle);
    matrix.at<double>(0, 1) = 0;
    matrix.at<double>(0, 2) = cos(angle);

    matrix.at<double>(1, 0) = 0;
    matrix.at<double>(1, 1) = 0;
    matrix.at<double>(1, 2) = 0;

    matrix.at<double>(2, 0) = -cos(angle);
    matrix.at<double>(2, 1) = 0;
    matrix.at<double>(2, 2) = -sin(angle);
    return matrix;
}

cv::Mat LeastSquaresPoly::DRz(double angle)
{
    cv::Mat matrix = cv::Mat::zeros(3, 3, CV_64F);

    matrix.at<double>(0, 0) = -sin(angle);
    matrix.at<double>(0, 1) = -cos(angle);
    matrix.at<double>(0, 2) = 0;

    matrix.at<double>(1, 0) = cos(angle);
    matrix.at<double>(1, 1) = -sin(angle);
    matrix.at<double>(1, 2) = 0;

    matrix.at<double>(2, 0) = 0;
    matrix.at<double>(2, 1) = 0;
    matrix.at<double>(2, 2) = 0;
    return matrix;
}

cv::Mat LeastSquaresPoly::CreatePoint(double x, double y, double z)
{
    cv::Mat trans(3, 1, CV_64F);
    trans.at<double>(0, 0) = x;
    trans.at<double>(1, 0) = y;
    trans.at<double>(2, 0) = z;
    return trans;
}

cv::Mat LeastSquaresPoly::CreatePoint(cv::Point3d Point)
{
    cv::Mat trans(3, 1, CV_64F);
    trans.at<double>(0, 0) = Point.x;
    trans.at<double>(1, 0) = Point.y;
    trans.at<double>(2, 0) = Point.z;
    return trans;
}

double LeastSquaresPoly::DF_Rx(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation, const cv::Mat& Point)
{
    double norm1 = getDistance(angleX, angleY, angleZ, pTranslation, Point);
    cv::Mat result2 = (DRx(angleX) * Ry(angleY) * Rz(angleZ)) * Point;
    cv::Point3d dPoint = cv::Point3d(result2);
    double norm2 = sqrt(dPoint.dot(dPoint));
    return norm1 * norm2;
}

double LeastSquaresPoly::DF_Ry(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation, const cv::Mat& Point)
{
    double norm1 = getDistance(angleX, angleY, angleZ, pTranslation, Point);
    cv::Mat result2 = (Rx(angleX) * DRy(angleY) * Rz(angleZ)) * Point;
    cv::Point3d dPoint = cv::Point3d(result2);
    double norm2 = sqrt(dPoint.dot(dPoint));
    return norm1 * norm2;
}

double LeastSquaresPoly::DF_Rz(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation, const cv::Mat& Point)
{
    double norm1 = getDistance(angleX, angleY, angleZ, pTranslation, Point);
    cv::Mat result2 = (Rx(angleX) * Ry(angleY) * DRz(angleZ)) * Point;
    cv::Point3d dPoint = cv::Point3d(result2);
    double norm2 = sqrt(dPoint.dot(dPoint));
    return norm1 * norm2;
}

double LeastSquaresPoly::DF_Translation(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation, const cv::Mat& Point)
{
    double norm1 = getDistance(angleX, angleY, angleZ, pTranslation, Point);
    return norm1;
}

double LeastSquaresPoly::getDistance(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation, const cv::Mat& Point)
{
    cv::Mat pointMat = (Rx(angleX) * Ry(angleY) * Rz(angleZ)) * Point + pTranslation;
    cv::Point3d tPoint = cv::Point3d(pointMat);
    double myClosest[3];
    double myPoint[3] = { tPoint.x, tPoint.y, tPoint.z };
    double signedDistance = implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(myPoint, myClosest);
    //return abs(signedDistance);

    return sqrt(tPoint.dot(tPoint));
}

cv::Point3d LeastSquaresPoly::ClosestPoint(const vtkSmartPointer<vtkPolyData>& surface, double point[3])
{
    double myClosest[3];
    double signedDistance = implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(point, myClosest);
    cv::Point3d closest = cv::Point3d(myClosest[0], myClosest[1], myClosest[2]);
    return closest;
}

double LeastSquaresPoly::computeFunction(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation)
{
    int tSize = source.size();
    double error = 0;
    for (int i = 0; i < tSize; i++)
    {
        cv::Mat Point = CreatePoint(source[i]);
        error += pow(getDistance(angleX, angleY, angleZ, pTranslation, Point), 2);
    }

    error = (error * 0.5) / double(tSize);
    return error;
}

cv::Mat LeastSquaresPoly::getGradientDiff(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation)
{
    double epsilon = 0.0001;
    cv::Point3d pointX = cv::Point3d(epsilon, 0, 0);
    cv::Point3d pointY = cv::Point3d(0, epsilon, 0);
    cv::Point3d pointZ = cv::Point3d(0, 0, epsilon);

    cv::Mat myTranslationX = pTranslation + CreatePoint(pointX);
    cv::Mat myTranslationY = pTranslation + CreatePoint(pointY);
    cv::Mat myTranslationZ = pTranslation + CreatePoint(pointZ);

    double myAngleX = angleX + epsilon;
    double myAngleY = angleY + epsilon;
    double myAngleZ = angleZ + epsilon;

    double aX = (computeFunction(myAngleX, angleY, angleZ, pTranslation) - computeFunction(angleX, angleY, angleZ, pTranslation)) / epsilon;
    double aY = (computeFunction(angleX, myAngleY, angleZ, pTranslation) - computeFunction(angleX, angleY, angleZ, pTranslation)) / epsilon;
    double aZ = (computeFunction(angleX, angleY, myAngleZ, pTranslation) - computeFunction(angleX, angleY, angleZ, pTranslation)) / epsilon;

    double X = (computeFunction(angleX, angleY, angleZ, myTranslationX) - computeFunction(angleX, angleY, angleZ, pTranslation)) / epsilon;
    double Y = (computeFunction(angleX, angleY, angleZ, myTranslationY) - computeFunction(angleX, angleY, angleZ, pTranslation)) / epsilon;
    double Z = (computeFunction(angleX, angleY, angleZ, myTranslationZ) - computeFunction(angleX, angleY, angleZ, pTranslation)) / epsilon;
    
    cv::Mat gradient(6, 1, CV_64F);

    gradient.at<double>(0, 0) = X;
    gradient.at<double>(1, 0) = Y;
    gradient.at<double>(2, 0) = Z;
    gradient.at<double>(3, 0) = aX;
    gradient.at<double>(4, 0) = aY;
    gradient.at<double>(5, 0) = aZ;

    return gradient;
}

cv::Mat LeastSquaresPoly::getGradient(double angleX, double angleY, double angleZ, const cv::Mat& pTranslation)
{
    int tSize = source.size();

    double translation = 0, aX = 0, aY = 0, aZ = 0;

    for (int i = 0; i < tSize; i++)
    {
        cv::Mat Point = CreatePoint(source[i]);
        translation += DF_Translation(angleX, angleY, angleZ, pTranslation, Point);
        aX += DF_Rx(angleX, angleY, angleZ, pTranslation, Point);
        aY += DF_Ry(angleX, angleY, angleZ, pTranslation, Point);
        aZ += DF_Rz(angleX, angleY, angleZ, pTranslation, Point);
    }

    cv::Mat gradient(6, 1, CV_64F);

    gradient.at<double>(0, 0) = translation;
    gradient.at<double>(1, 0) = translation;
    gradient.at<double>(2, 0) = translation;
    gradient.at<double>(3, 0) = aX;
    gradient.at<double>(4, 0) = aY;
    gradient.at<double>(5, 0) = aZ;

    return (gradient / double(tSize));
}

double LeastSquaresPoly::LeastSquares(cv::Mat& data, double learningRate, int iterations)
{
    int tSize = source.size();
    cv::Mat translation(3, 1, CV_64F);
    double error = 0;
    for (int i = 0; i < iterations; i++)
    {
        translation.at<double>(0, 0) = data.at<double>(0, 0);
        translation.at<double>(1, 0) = data.at<double>(1, 0);
        translation.at<double>(2, 0) = data.at<double>(2, 0);

        double angleX = data.at<double>(3, 0);
        double angleY = data.at<double>(4, 0);
        double angleZ = data.at<double>(5, 0);
        
        cv::Mat gradient = getGradient(angleX, angleY, angleZ, translation);
        cv::Mat gradientDiff = getGradientDiff(angleX, angleY, angleZ, translation);

        std::cout << gradient << std::endl;
        std::cout << gradientDiff << std::endl;
        //std::cout << computeFunction(angleX, angleY, angleZ, translation) << std::endl;
        std::cout << "*****************************************************************" << std::endl;

        data = data - learningRate * gradientDiff;
    }
    return error;
}