#include "LeastSquaresICP.hpp"

///////////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>
#include <vtkNamedColors.h>
#include <vtkPolyDataMapper.h>
#include <algorithm>
#include <random>

using namespace SPINE::SEGMENTATION;

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

////////////////////////////////////////////////////////////////////////////////////


LeastSquaresICP::LeastSquaresICP(const std::vector<PointTypeITK>& sourcePoints)
{
    for (int i = 0; i < sourcePoints.size(); i++)
    {
        source.push_back(cv::Point3d(sourcePoints[i][0], sourcePoints[i][1], sourcePoints[i][2]));
    }

    chi2 = 0.75;
    maxError = 0.3;
}

LeastSquaresICP::LeastSquaresICP(const std::vector<cv::Point3d>& sourcePoints)
{
    for (int i = 0; i < sourcePoints.size(); i++)
    {
        source.push_back(sourcePoints[i]);
    }

    chi2 = 0.75;
    maxError = 0.3;
}

void LeastSquaresICP::setChi2(double pChi2)
{
    chi2 = pChi2;
}

void LeastSquaresICP::setMaxError(double pMaxError)
{
    maxError = pMaxError;
}

double LeastSquaresICP::getChi2() const
{
    return chi2;
}

double LeastSquaresICP::getMaxError() const
{
    return maxError;
}

void LeastSquaresICP::shuffleSource()
{
    auto rng = std::default_random_engine{};
    std::shuffle(source.begin(), source.end(), rng);
}

cv::Mat LeastSquaresICP::Rx(double angle)
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

cv::Mat LeastSquaresICP::Ry(double angle)
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

cv::Mat LeastSquaresICP::Rz(double angle)
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

cv::Mat LeastSquaresICP::DRx(double angle)
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

cv::Mat LeastSquaresICP::DRy(double angle)
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

cv::Mat LeastSquaresICP::DRz(double angle)
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

cv::Mat LeastSquaresICP::CreatePoint(double x, double y, double z)
{
    cv::Mat trans(3, 1, CV_64F);
    trans.at<double>(0, 0) = x;
    trans.at<double>(1, 0) = y;
    trans.at<double>(2, 0) = z;
    return trans;
}

cv::Mat LeastSquaresICP::CreatePoint(cv::Point3d Point)
{
    cv::Mat trans(3, 1, CV_64F);
    trans.at<double>(0, 0) = Point.x;
    trans.at<double>(1, 0) = Point.y;
    trans.at<double>(2, 0) = Point.z;
    return trans;
}

cv::Mat LeastSquaresICP::DF_Rx(double angleX, double angleY, double angleZ, const cv::Mat& Point)
{
    cv::Mat result = (DRx(angleX) * Ry(angleY) * Rz(angleZ)) * Point;
    return result;
}

cv::Mat LeastSquaresICP::DF_Ry(double angleX, double angleY, double angleZ, const cv::Mat& Point)
{
    cv::Mat result = (Rx(angleX) * DRy(angleY) * Rz(angleZ)) * Point;
    return result;
}

cv::Mat LeastSquaresICP::DF_Rz(double angleX, double angleY, double angleZ, const cv::Mat& Point)
{
    cv::Mat result = (Rx(angleX) * Ry(angleY) * DRz(angleZ)) * Point;
    return result;
}

cv::Mat LeastSquaresICP::Jacobian(const cv::Mat& data, const cv::Mat& Point)
{
    double angleX = data.at<double>(3, 0);
    double angleY = data.at<double>(4, 0);
    double angleZ = data.at<double>(5, 0);

    cv::Mat H1, H2, J;
    cv::Mat translation = cv::Mat::eye(3, 3, CV_64F);
    cv::Mat alphaX = DF_Rx(angleX, angleY, angleZ, Point);
    cv::Mat alphaY = DF_Ry(angleX, angleY, angleZ, Point);
    cv::Mat alphaZ = DF_Rz(angleX, angleY, angleZ, Point);

    cv::hconcat(translation, alphaX, H1);
    cv::hconcat(H1, alphaY, H2);
    cv::hconcat(H2, alphaZ, J);

    return J;
}

cv::Mat LeastSquaresICP::SquareError(const cv::Mat& data, const cv::Mat& source, const cv::Mat& target)
{
    cv::Mat translation(3, 1, CV_64F);
    translation.at<double>(0, 0) = data.at<double>(0, 0);
    translation.at<double>(1, 0) = data.at<double>(1, 0);
    translation.at<double>(2, 0) = data.at<double>(2, 0);

    cv::Mat result = (Rx(data.at<double>(3, 0)) * Ry(data.at<double>(4, 0)) * Rz(data.at<double>(5, 0))) * source + translation;
    cv::Mat error = result - target;
    return error;
}

LeastSquaresICP::GaussNewton LeastSquaresICP::GetSystem(const std::vector<PointTypeITK>& target, const cv::Mat& data, int posBegin, int posEnd, double lambda)
{
    cv::Mat A = cv::Mat::zeros(6, 6, CV_64F);
    cv::Mat B = cv::Mat::zeros(6, 1, CV_64F);
    double chi = 0;
    double squareError;
    double localError = 0;
    int tSize = source.size();

    for (int i = 0; i < tSize; i++)
    {
        cv::Mat sourcePoint = CreatePoint(source[i]);
        cv::Mat targetPoint = CreatePoint(target[i][0], target[i][1], target[i][2]);
        cv::Mat error = SquareError(data, sourcePoint, targetPoint);
		squareError = error.dot(error);

		if (i >= posBegin && i < posEnd)
		{
			cv::Mat Jac = Jacobian(data, sourcePoint);

			cv::Mat squareJac = (Jac.t())*Jac;
			cv::Mat d0 = squareJac.diag(0);
			cv::Mat diagonal = cv::Mat::diag(d0);

			A = A + (squareJac + lambda * diagonal);
			B = B + (Jac.t())*error;

			localError += squareError;
		}

        chi = chi + squareError;
    }
    return GaussNewton(A, -B, chi, localError);
}

LeastSquaresICP::GaussNewton LeastSquaresICP::GetSystem(const std::vector<cv::Point3d>& target, const cv::Mat& data, int posBegin, int posEnd, double lambda)
{
    cv::Mat A = cv::Mat::zeros(6, 6, CV_64F);
    cv::Mat B = cv::Mat::zeros(6, 1, CV_64F);
    double chi = 0;
    double squareError;
    double localError = 0;
    int tSize = source.size();

    for (int i = 0; i < tSize; i++)
    {
        cv::Mat sourcePoint = CreatePoint(source[i]);
        cv::Mat targetPoint = CreatePoint(target[i]);
        cv::Mat error = SquareError(data, sourcePoint, targetPoint);
		squareError = error.dot(error);

		if ( i >= posBegin && i < posEnd )
		{
			cv::Mat Jac = Jacobian(data, sourcePoint);

			cv::Mat squareJac = (Jac.t())*Jac;
			cv::Mat d0 = squareJac.diag(0);
			cv::Mat diagonal = cv::Mat::diag(d0);

			A = A + (squareJac + lambda * diagonal);
			B = B + (Jac.t())*error;

			localError += squareError;
		}

        chi = chi + squareError;
    }

    return GaussNewton(A, -B, chi, localError);
}

std::vector<cv::Point3d> LeastSquaresICP::GetCorrespondence(const cv::Mat& target, const cv::Mat& data)
{
    std::vector<cv::Point3d> result;

    double angleX = data.at<double>(3, 0);
    double angleY = data.at<double>(4, 0);
    double angleZ = data.at<double>(5, 0);

    cv::Mat translation(3, 1, CV_64F);
    translation.at<double>(0, 0) = data.at<double>(0, 0);
    translation.at<double>(1, 0) = data.at<double>(1, 0);
    translation.at<double>(2, 0) = data.at<double>(2, 0);

	cv::Mat sourceMat(source.size(), 3, CV_64F);

    for (int i = 0; i < source.size(); i++)
    {
        cv::Mat pointMat = CreatePoint(source[i]);
        cv::Mat newPointMat = (Rx(angleX) * Ry(angleY) * Rz(angleZ)) * pointMat + translation;
        cv::Point3d newPoint = cv::Point3d(newPointMat);

		sourceMat.at<double>(i, 0) = newPoint.x;
		sourceMat.at<double>(i, 1) = newPoint.y;
		sourceMat.at<double>(i, 2) = newPoint.z;
    }

	cv::FlannBasedMatcher matcher(new cv::flann::KDTreeIndexParams(1));
	std::vector<std::vector<cv::DMatch>> knnMatches;
	matcher.knnMatch(sourceMat, target, knnMatches, 1); // 1 KNN

	for (size_t i = 0; i < knnMatches.size(); ++i) {
		int closestIdx = knnMatches[i][0].trainIdx;
		//float closestDist = knnMatches[i][0].distance;
		cv::Point3d tempPoint = cv::Point3d(target.at<double>(closestIdx, 0), target.at<double>(closestIdx, 1), target.at<double>(closestIdx, 2));
		result.push_back(tempPoint);
	}


    return result;
}

double LeastSquaresICP::LeastSquares(const cv::Mat& targetOnCT, cv::Mat& data, int iterations)
{
    std::vector<cv::Point3d> target = GetCorrespondence(targetOnCT, data);
    double angleX, angleY, angleZ;

    bool finish = false;

    double lambda = 0.01;
    double currentError, beforeError = -1, totalError;
    double maxLambda = 1000.0;

    int batch = 3;
    int step = source.size() / batch;

    cv::Mat dataTemp(6, 1, CV_64F);
    double bestError = -1;
	double tSize = source.size();

    for (int i = 0; i < iterations && finish == false; i++)
    {

        for (int j = 0; j < batch; j++)
        {
            int posA, posB;
            posA = step * j;
            posB = posA + step;

            if (j == batch - 1)
            {
                posB = source.size();
            }

            GaussNewton resultInfo = GetSystem(target, data, posA, posB, lambda);
            currentError = resultInfo.localError;
			totalError = sqrt(resultInfo.totalError / tSize);

			if (bestError < 0 || totalError < bestError)
			{
				bestError = totalError;

				dataTemp.at<double>(0, 0) = data.at<double>(0, 0);
				dataTemp.at<double>(1, 0) = data.at<double>(1, 0);
				dataTemp.at<double>(2, 0) = data.at<double>(2, 0);
				dataTemp.at<double>(3, 0) = data.at<double>(3, 0);
				dataTemp.at<double>(4, 0) = data.at<double>(4, 0);
				dataTemp.at<double>(5, 0) = data.at<double>(5, 0);
			}

			if (bestError < maxError)
			{
				finish = true;
				break;
			}

            if (beforeError < 0)
            {
                beforeError = currentError;
            }
            else
            {
                if (currentError < beforeError)
                {
                    lambda = lambda / 3.0;
                }
                else if (currentError > beforeError)
                {
                    lambda = lambda * 2.0;
                    if (lambda > maxLambda)
                    {
                        lambda = maxLambda;
                    }
                }
                beforeError = currentError;
            }

            cv::Mat dx;

            bool result = cv::solve(resultInfo.A, resultInfo.B, dx);

            data = data + dx;

            angleX = data.at<double>(3, 0);
            angleY = data.at<double>(4, 0);
            angleZ = data.at<double>(5, 0);

            data.at<double>(3, 0) = atan2(sin(angleX), cos(angleX));
            data.at<double>(4, 0) = atan2(sin(angleY), cos(angleY));
            data.at<double>(5, 0) = atan2(sin(angleZ), cos(angleZ));

            /*if (resultInfo.totalError < chi2 && resultInfo.localError < maxError)
            {
                finish = true;
                break;
            }*/

        }

		/*std::sort(tempResult.begin(), tempResult.end(), [](const BatchResult &a, const BatchResult &b) { return (a.error < b.error); });
		int pos = batch - 2;
		if (bestError < 0 || tempResult[pos].error < bestError)
		{
			bestError = tempResult[pos].error;

			dataTemp.at<double>(0, 0) = tempResult[pos].data.at<double>(0, 0);
			dataTemp.at<double>(1, 0) = tempResult[pos].data.at<double>(1, 0);
			dataTemp.at<double>(2, 0) = tempResult[pos].data.at<double>(2, 0);
			dataTemp.at<double>(3, 0) = tempResult[pos].data.at<double>(3, 0);
			dataTemp.at<double>(4, 0) = tempResult[pos].data.at<double>(4, 0);
			dataTemp.at<double>(5, 0) = tempResult[pos].data.at<double>(5, 0);
		}

		tempResult.clear();*/

		//if (bestError < maxError)
		//{
		//	finish = true;
		//}

		if (finish == true)
		{
			continue;
		}

        shuffleSource();
        target.clear();
        target = GetCorrespondence(targetOnCT, data);
    }

    /* cv::Mat translation(3, 1, CV_64F);
     translation.at<double>(0, 0) = data.at<double>(0, 0);
     translation.at<double>(1, 0) = data.at<double>(1, 0);
     translation.at<double>(2, 0) = data.at<double>(2, 0);

     for (int i = 0; i < source.size(); i++)
     {
         cv::Mat transformPoint(3, 1, CV_64F);

         transformPoint.at<double>(0, 0) = source[i].x;
         transformPoint.at<double>(1, 0) = source[i].y;
         transformPoint.at<double>(2, 0) = source[i].z;

         cv::Mat result = (Rx(data.at<double>(3, 0)) * Ry(data.at<double>(4, 0)) * Rz(data.at<double>(5, 0))) * transformPoint + translation;

         cv::Point3d myLastPoint = cv::Point3d(result);

         double pnt[3];

         pnt[0] = myLastPoint.x;
         pnt[1] = myLastPoint.y;
         pnt[2] = myLastPoint.z;

         ClosestPoint(surface, pnt);
     }*/

    data.at<double>(0, 0) = dataTemp.at<double>(0, 0);
    data.at<double>(1, 0) = dataTemp.at<double>(1, 0);
    data.at<double>(2, 0) = dataTemp.at<double>(2, 0);
    data.at<double>(3, 0) = dataTemp.at<double>(3, 0);
    data.at<double>(4, 0) = dataTemp.at<double>(4, 0);
    data.at<double>(5, 0) = dataTemp.at<double>(5, 0);

    return bestError;
}

double LeastSquaresICP::LeastSquaresRandomInit(const cv::Mat& targetOnCT, cv::Mat& data, int iterations)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distr_translation(-1, 1);
    std::uniform_real_distribution<> distr_rotation(-0.06, 0.06);

    cv::Mat dataTemp(6, 1, CV_64F);
    cv::Mat dataInit(6, 1, CV_64F);
    cv::Mat randomValues(6, 1, CV_64F);

    dataInit.at<double>(0, 0) = data.at<double>(0, 0);
    dataInit.at<double>(1, 0) = data.at<double>(1, 0);
    dataInit.at<double>(2, 0) = data.at<double>(2, 0);
    dataInit.at<double>(3, 0) = data.at<double>(3, 0);
    dataInit.at<double>(4, 0) = data.at<double>(4, 0);
    dataInit.at<double>(5, 0) = data.at<double>(5, 0);

    double error = LeastSquares(targetOnCT, data, iterations);

    for (int i = 0; i < 10; i++)
    {

        dataTemp.at<double>(0, 0) = dataInit.at<double>(0, 0);
        dataTemp.at<double>(1, 0) = dataInit.at<double>(1, 0);
        dataTemp.at<double>(2, 0) = dataInit.at<double>(2, 0);
        dataTemp.at<double>(3, 0) = dataInit.at<double>(3, 0);
        dataTemp.at<double>(4, 0) = dataInit.at<double>(4, 0);
        dataTemp.at<double>(5, 0) = dataInit.at<double>(5, 0);

        randomValues.at<double>(0, 0) = distr_translation(gen);
        randomValues.at<double>(1, 0) = distr_translation(gen);
        randomValues.at<double>(2, 0) = distr_translation(gen);

        randomValues.at<double>(3, 0) = distr_rotation(gen);
        randomValues.at<double>(4, 0) = distr_rotation(gen);
        randomValues.at<double>(5, 0) = distr_rotation(gen);

        dataTemp += randomValues;

        double errorTemp = LeastSquares(targetOnCT, dataTemp, iterations);

        if (errorTemp < error)
        {
            error = errorTemp;

            data.at<double>(0, 0) = dataTemp.at<double>(0, 0);
            data.at<double>(1, 0) = dataTemp.at<double>(1, 0);
            data.at<double>(2, 0) = dataTemp.at<double>(2, 0);
            data.at<double>(3, 0) = dataTemp.at<double>(3, 0);
            data.at<double>(4, 0) = dataTemp.at<double>(4, 0);
            data.at<double>(5, 0) = dataTemp.at<double>(5, 0);
        }

    }

    return error;
}

cv::Mat LeastSquaresICP::GetRotationAnglesXYZ(const std::vector<cv::Point3d>& threeVectorsSource, const std::vector<cv::Point3d>& threeVectorstarget, cv::Mat& data)
{
    std::vector<cv::Point3d> sources = threeVectorsSource;
    std::vector<cv::Point3d> targets = threeVectorstarget;

    cv::Mat sourceMatrix = cv::Mat(sources.size(), 3, CV_64F, sources.data());
    cv::Mat targetMatrix = cv::Mat(targets.size(), 3, CV_64F, targets.data());

    cv::Mat inverse = (sourceMatrix.t()).inv();
    cv::Mat rotationCV = (targetMatrix.t()) * inverse;

    Eigen::Matrix3d rotationEIGEN(3, 3);

    rotationEIGEN(0, 0) = rotationCV.at<double>(0, 0);
    rotationEIGEN(1, 0) = rotationCV.at<double>(1, 0);
    rotationEIGEN(2, 0) = rotationCV.at<double>(2, 0);

    rotationEIGEN(0, 1) = rotationCV.at<double>(0, 1);
    rotationEIGEN(1, 1) = rotationCV.at<double>(1, 1);
    rotationEIGEN(2, 1) = rotationCV.at<double>(2, 1);

    rotationEIGEN(0, 2) = rotationCV.at<double>(0, 2);
    rotationEIGEN(1, 2) = rotationCV.at<double>(1, 2);
    rotationEIGEN(2, 2) = rotationCV.at<double>(2, 2);

    Eigen::Vector3d rot = rotationEIGEN.eulerAngles(0, 1, 2);

    data.at<double>(0, 0) = rot(0);
    data.at<double>(1, 0) = rot(1);
    data.at<double>(2, 0) = rot(2);

    return rotationCV;
}

cv::Mat LeastSquaresICP::GetRotationMatrix(double angleX, double angleY, double angleZ)
{
    return (Rx(angleX) * Ry(angleY) * Rz(angleZ));
}
