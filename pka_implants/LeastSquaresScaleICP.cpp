#include "LeastSquaresScaleICP.hpp"
#include <Eigen/Dense>

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
#include <random>

using namespace PKA::IMPLANTS;

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

void showProcess(vtkSmartPointer<vtkPolyData> poly, std::vector<cv::Point3d> pPoints)
{
    vtkSmartPointer<vtkPolyData> surfaceJoin = MergePolyWithSphere(poly, pPoints);
    ShowPolyData(surfaceJoin, poly);
}

////////////////////////////////////////////////////////////////////////////////////


LeastSquaresScaleICP::LeastSquaresScaleICP(const std::vector<cv::Point3d>& sourcePoints)
{
    aveSource = cv::Point3d(0, 0, 0);
    for (int i = 0; i < sourcePoints.size(); i++)
    {
        source.push_back(cv::Point3d(sourcePoints[i]));
        aveSource = aveSource + source[i];
    }

    if (source.size() > 0)
    {
        aveSource = aveSource / double(source.size());

        for (int i = 0; i < source.size(); i++)
        {
            cv::Point3d temp = source[i] - aveSource;
            centerSource.push_back(temp);
        }
    }

    chi2 = 0.75;
    maxError = 1.0;
}

void LeastSquaresScaleICP::shuffleCenterSource()
{
    auto rng = std::default_random_engine{};
    std::shuffle(centerSource.begin(), centerSource.end(), rng);
}

void LeastSquaresScaleICP::setChi2(double pChi2)
{
    chi2 = pChi2;
}

void LeastSquaresScaleICP::setMaxError(double pMaxError)
{
    maxError = pMaxError;
}

double LeastSquaresScaleICP::getChi2() const
{
    return chi2;
}

double LeastSquaresScaleICP::getMaxError() const
{
    return maxError;
}

cv::Mat LeastSquaresScaleICP::Rx(double angle)
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

cv::Mat LeastSquaresScaleICP::Ry(double angle)
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

cv::Mat LeastSquaresScaleICP::Rz(double angle)
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

cv::Mat LeastSquaresScaleICP::DRx(double angle)
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

cv::Mat LeastSquaresScaleICP::DRy(double angle)
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

cv::Mat LeastSquaresScaleICP::DRz(double angle)
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

cv::Mat LeastSquaresScaleICP::CreatePoint(double x, double y, double z)
{
    cv::Mat trans(3, 1, CV_64F);
    trans.at<double>(0, 0) = x;
    trans.at<double>(1, 0) = y;
    trans.at<double>(2, 0) = z;
    return trans;
}

cv::Mat LeastSquaresScaleICP::CreatePoint(cv::Point3d Point)
{
    cv::Mat trans(3, 1, CV_64F);
    trans.at<double>(0, 0) = Point.x;
    trans.at<double>(1, 0) = Point.y;
    trans.at<double>(2, 0) = Point.z;
    return trans;
}

cv::Mat LeastSquaresScaleICP::DF_Rx(double angleX, double angleY, double angleZ, const cv::Mat& Point)
{
    cv::Mat result = (DRx(angleX) * Ry(angleY) * Rz(angleZ)) * Point;
    return result;
}

cv::Mat LeastSquaresScaleICP::DF_Ry(double angleX, double angleY, double angleZ, const cv::Mat& Point)
{
    cv::Mat result = (Rx(angleX) * DRy(angleY) * Rz(angleZ)) * Point;
    return result;
}

cv::Mat LeastSquaresScaleICP::DF_Rz(double angleX, double angleY, double angleZ, const cv::Mat& Point)
{
    cv::Mat result = (Rx(angleX) * Ry(angleY) * DRz(angleZ)) * Point;
    return result;
}

cv::Mat LeastSquaresScaleICP::JacobianScale(const cv::Mat& data, const cv::Mat& Point)
{
    double angleX = data.at<double>(3, 0);
    double angleY = data.at<double>(4, 0);
    double angleZ = data.at<double>(5, 0);

    double scale = data.at<double>(6, 0);

    cv::Mat H1, H2, J;
    cv::Mat translationDf = (cv::Mat::eye(3, 3, CV_64F));

    cv::Mat alphaX = scale * DF_Rx(angleX, angleY, angleZ, Point);
    cv::Mat alphaY = scale * DF_Ry(angleX, angleY, angleZ, Point);
    cv::Mat alphaZ = scale * DF_Rz(angleX, angleY, angleZ, Point);

    //cv::Mat scaleDf = -(1.0 / (scale * scale)) * ((Rx(angleX) * Ry(angleY) * Rz(angleZ)) * Point) - (2.0 / (scale * scale * scale)) * translation + (2.0 / (scale * scale * scale)) * PointBi;

    cv::hconcat(translationDf, alphaX, H1);
    cv::hconcat(H1, alphaY, H2);
    cv::hconcat(H2, alphaZ, J);

    return J;
}


cv::Mat LeastSquaresScaleICP::SquareErrorScale(const cv::Mat& data, const cv::Mat& source, const cv::Mat& target)
{
    cv::Mat translation(3, 1, CV_64F);

    translation.at<double>(0, 0) = data.at<double>(0, 0);
    translation.at<double>(1, 0) = data.at<double>(1, 0);
    translation.at<double>(2, 0) = data.at<double>(2, 0);

    double scale = data.at<double>(6, 0);

    cv::Mat rotation = Rx(data.at<double>(3, 0)) * Ry(data.at<double>(4, 0)) * Rz(data.at<double>(5, 0));

    cv::Mat result = scale * (rotation * source) + translation;
    cv::Mat error = result - target;
    return error;
}

void LeastSquaresScaleICP::GetScale(const std::vector<cv::Point3d>& target, cv::Mat& data)
{
    cv::Mat rotation;
    rotation = Rx(data.at<double>(3, 0)) * Ry(data.at<double>(4, 0)) * Rz(data.at<double>(5, 0));

    cv::Mat translation(3, 1, CV_64F);
    translation.at<double>(0, 0) = data.at<double>(0, 0);
    translation.at<double>(1, 0) = data.at<double>(1, 0);
    translation.at<double>(2, 0) = data.at<double>(2, 0);

    double scale = data.at<double>(6, 0);

    int tSize = centerSource.size();

    double add1 = 0;
    double add2 = 0;

    for (int i = 0; i < tSize; i++)
    {
        cv::Mat movePointRotation = scale * (rotation * CreatePoint(centerSource[i]));
        cv::Mat movePointTarget = translation - CreatePoint(target[i]);

        cv::Mat varA = movePointRotation.t() * movePointTarget;
        cv::Mat varB = movePointTarget.t() * movePointTarget;

        add1 = add1 + varA.at<double>(0, 0);
        add2 = add2 + varB.at<double>(0, 0);
    }

    data.at<double>(6, 0) = scale * (-(add2 / add1));
}

LeastSquaresScaleICP::GaussNewton LeastSquaresScaleICP::GetSystemScale(const std::vector<cv::Point3d>& target, const cv::Mat& data, int posBegin, int posEnd, double lambda)
{
    cv::Mat A = cv::Mat::zeros(6, 6, CV_64F);
    cv::Mat B = cv::Mat::zeros(6, 1, CV_64F);
    double chi = 0;
    double squareError;
    double localError = 0;

    for (int i = posBegin; i < posEnd; i++)
    {
        cv::Mat sourcePoint = CreatePoint(centerSource[i]);
        cv::Mat targetPoint = CreatePoint(target[i]);
        cv::Mat error = SquareErrorScale(data, sourcePoint, targetPoint);
		squareError = error.dot(error);

		if (i >= posBegin && i < posEnd)
		{
			cv::Mat Jac = JacobianScale(data, sourcePoint);

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

cv::Point3d LeastSquaresScaleICP::ClosestPoint(const vtkSmartPointer<vtkPolyData>& surface, double point[3])
{
    vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
    implicitPolyDataDistance->SetInput(surface);
    double myClosest[3];
    double signedDistance = implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(point, myClosest);
    cv::Point3d closest = cv::Point3d(myClosest[0], myClosest[1], myClosest[2]);
    return closest;
}

std::vector<cv::Point3d> LeastSquaresScaleICP::GetCorrespondenceScale(const vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance, const cv::Mat& data)
{
    std::vector<cv::Point3d> result;

    double angleX = data.at<double>(3, 0);
    double angleY = data.at<double>(4, 0);
    double angleZ = data.at<double>(5, 0);

    cv::Mat translation(3, 1, CV_64F);
    translation.at<double>(0, 0) = data.at<double>(0, 0);
    translation.at<double>(1, 0) = data.at<double>(1, 0);
    translation.at<double>(2, 0) = data.at<double>(2, 0);

    double scale = data.at<double>(6, 0);

    for (int i = 0; i < centerSource.size(); i++)
    {
        cv::Mat pointMat = CreatePoint(centerSource[i]);
        cv::Mat newPointMat = scale * ((Rx(angleX) * Ry(angleY) * Rz(angleZ)) * pointMat) + translation;
        cv::Point3d newPoint = cv::Point3d(newPointMat);

        double pnt[3];
        pnt[0] = newPoint.x;
        pnt[1] = newPoint.y;
        pnt[2] = newPoint.z;

        double myClosest[3];
        implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);

        result.push_back(cv::Point3d(myClosest[0], myClosest[1], myClosest[2]));
    }

    return result;
}

cv::Mat LeastSquaresScaleICP::GetRotationAnglesXYZ(const std::vector<cv::Point3d>& threeVectorsSource, const std::vector<cv::Point3d>& threeVectorstarget, cv::Mat& data)
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

cv::Mat LeastSquaresScaleICP::GetRotationMatrix(double angleX, double angleY, double angleZ)
{
    return (Rx(angleX) * Ry(angleY) * Rz(angleZ));
}

double LeastSquaresScaleICP::LeastSquaresScale(const vtkSmartPointer<vtkPolyData>& surface, cv::Mat& data, int iterations)
{
    vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
    implicitPolyDataDistance->SetInput(surface);

    cv::Mat myRotation = Rx(data.at<double>(3, 0)) * Ry(data.at<double>(4, 0)) * Rz(data.at<double>(5, 0));
    cv::Mat myRest = (myRotation * CreatePoint(aveSource));

    data.at<double>(0, 0) = data.at<double>(0, 0) + myRest.at<double>(0, 0);
    data.at<double>(1, 0) = data.at<double>(1, 0) + myRest.at<double>(1, 0);
    data.at<double>(2, 0) = data.at<double>(2, 0) + myRest.at<double>(2, 0);

    std::vector<cv::Point3d> target = GetCorrespondenceScale(implicitPolyDataDistance, data);
    double angleX, angleY, angleZ;

    bool finish = false;
    double lambda = 0.01;
    double maxLambda = 1000.0;
    double currentError, beforeError = -1, totalError;

    int batch = 3;
    int step = source.size() / batch;

    cv::Mat dataTemp(7, 1, CV_64F);
    double bestError = -1;
	double tSize = source.size();

    //////////////////////////////////////////////////////////
    
    /*cv::Mat translation1(3, 1, CV_64F);
    translation1.at<double>(0, 0) = data.at<double>(0, 0);
    translation1.at<double>(1, 0) = data.at<double>(1, 0);
    translation1.at<double>(2, 0) = data.at<double>(2, 0);

    cv::Mat finalRest1 = data.at<double>(6, 0) * ((Rx(data.at<double>(3, 0)) * Ry(data.at<double>(4, 0)) * Rz(data.at<double>(5, 0))) * CreatePoint(aveSource));
    translation1 = translation1 - finalRest1;
    std::vector<cv::Point3d> mySource1;
    for (int j = 0; j < source.size(); j++)
    {
        cv::Mat myPoint = CreatePoint(source[j]);
        cv::Mat eval = data.at<double>(6, 0) * ((Rx(data.at<double>(3, 0)) * Ry(data.at<double>(4, 0)) * Rz(data.at<double>(5, 0))) * myPoint) + translation1;
        mySource1.push_back(cv::Point3d(eval));
    }

    showProcess(surface, mySource1);*/

    ///////////////////////////////////////

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

            GaussNewton resultInfo = GetSystemScale(target, data, posA, posB, lambda);
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
				dataTemp.at<double>(6, 0) = data.at<double>(6, 0);
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

            data.at<double>(0, 0) = data.at<double>(0, 0) + dx.at<double>(0, 0);
            data.at<double>(1, 0) = data.at<double>(1, 0) + dx.at<double>(1, 0);
            data.at<double>(2, 0) = data.at<double>(2, 0) + dx.at<double>(2, 0);

            data.at<double>(3, 0) = data.at<double>(3, 0) + dx.at<double>(3, 0);
            data.at<double>(4, 0) = data.at<double>(4, 0) + dx.at<double>(4, 0);
            data.at<double>(5, 0) = data.at<double>(5, 0) + dx.at<double>(5, 0);

            angleX = data.at<double>(3, 0);
            angleY = data.at<double>(4, 0);
            angleZ = data.at<double>(5, 0);

            data.at<double>(3, 0) = atan2(sin(angleX), cos(angleX));
            data.at<double>(4, 0) = atan2(sin(angleY), cos(angleY));
            data.at<double>(5, 0) = atan2(sin(angleZ), cos(angleZ));

            target.clear();
            target = GetCorrespondenceScale(implicitPolyDataDistance, data);

            GetScale(target, data);
        }

        if (finish == true)
        {
            continue;
        }

        shuffleCenterSource();
        target.clear();
        target = GetCorrespondenceScale(implicitPolyDataDistance, data);
    }

    data.at<double>(0, 0) = dataTemp.at<double>(0, 0);
    data.at<double>(1, 0) = dataTemp.at<double>(1, 0);
    data.at<double>(2, 0) = dataTemp.at<double>(2, 0);
    data.at<double>(3, 0) = dataTemp.at<double>(3, 0);
    data.at<double>(4, 0) = dataTemp.at<double>(4, 0);
    data.at<double>(5, 0) = dataTemp.at<double>(5, 0);
    data.at<double>(6, 0) = dataTemp.at<double>(6, 0);

    double myScale = data.at<double>(6, 0);

    cv::Mat finalRest = myScale * ((Rx(data.at<double>(3, 0)) * Ry(data.at<double>(4, 0)) * Rz(data.at<double>(5, 0))) * CreatePoint(aveSource));

    data.at<double>(0, 0) = data.at<double>(0, 0) - finalRest.at<double>(0, 0);
    data.at<double>(1, 0) = data.at<double>(1, 0) - finalRest.at<double>(1, 0);
    data.at<double>(2, 0) = data.at<double>(2, 0) - finalRest.at<double>(2, 0);

    //////////////////////////////////////////////////////////

    
    /*cv::Mat translation(3, 1, CV_64F);
    translation.at<double>(0, 0) = data.at<double>(0, 0);
    translation.at<double>(1, 0) = data.at<double>(1, 0);
    translation.at<double>(2, 0) = data.at<double>(2, 0);

    std::vector<cv::Point3d> mySource;
    for (int j = 0; j < source.size(); j++)
    {
        cv::Mat myPoint = CreatePoint(source[j]);
        cv::Mat eval = myScale * ((Rx(data.at<double>(3, 0)) * Ry(data.at<double>(4, 0)) * Rz(data.at<double>(5, 0))) * myPoint) + translation;
        mySource.push_back(cv::Point3d(eval));
    }

    showProcess(surface, mySource);
*/
    ///////////////////////////////////////

    return bestError;
}
