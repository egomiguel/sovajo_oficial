
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
#include <iostream>
#include "FemurImplant.hpp"
#include "ImplantsException.hpp"

using namespace TKA::IMPLANTS;

FemurImplant::FemurImplant()
{
	isInit = false;
}

void FemurImplant::fixNormalVectorC()
{
	/*Point midPointCenter = (P1 + P2) / 2.0;
	Point midPoint1 = (midPointCenter + cortexPoint) / 2.0;
	Point normalVectorPlaneC = planeC.getNormalVector();
	Point basePoint = midPoint1 - normalVectorPlaneC;

	Line intereceptionLine(basePoint, midPoint1);

	Point iterceptionPoint = planeC.getInterceptionLinePoint(intereceptionLine);
	Point specialNormalVector = midPoint1 - iterceptionPoint;
	planeC.fixNormalVector(specialNormalVector);*/

    Point cortexProj = planeC.getProjectionPoint(cortexPoint);
    Point specialNormalVector = cortexPoint - cortexProj;
    planeC.fixNormalVector(specialNormalVector);
}

void FemurImplant::init(const Plane& A, const Plane& B, const Plane& C, const Plane& D, const Plane& E,
    const Point& P3, const Point& P4, const Point& cortexPoint, const vtkSmartPointer<vtkPolyData> implantModel,
    const std::vector<PointTypeITK>& kneeCapPath, const FemurImplantInfo& pImplantInfo)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_FEMUR_IMPLANT;
    }
    this->planeA = A;
    this->planeB = B;
    this->planeC = C;
    this->planeD = D;
    this->planeE = E;
    this->P3 = P3;
    this->P4 = P4;
    this->cortexPoint = cortexPoint;
    this->implantModel = implantModel;
    this->mImplantInfo = pImplantInfo;

    Point midPoint = (P3 + P4) / 2.0;
    Point directVector = P4 - P3;
    midPlane.init(directVector, midPoint);

    fixNormalVectorC();

    FixKneeCapPath(kneeCapPath);

    isInit = true;
}

FemurImplant::FemurImplant(const FemurImplant& pImplant)
{
    this->planeA = pImplant.planeA;
    this->planeB = pImplant.planeB;
    this->planeC = pImplant.planeC;
    this->planeD = pImplant.planeD;
    this->planeE = pImplant.planeE;
    this->P3 = pImplant.P3;
    this->P4 = pImplant.P4;
    this->cortexPoint = pImplant.cortexPoint;
    this->implantModel = pImplant.implantModel;
    this->isInit = pImplant.isInit;
    this->midPlane = pImplant.midPlane;
    this->mKneeCapPath = pImplant.mKneeCapPath;
    this->mImplantInfo = pImplant.mImplantInfo;
}

Point FemurImplant::getMidPlaneInterceptionPoint() const
{
    Line temp = Line::makeLineWithPoints(P3, P4);
    return midPlane.getInterceptionLinePoint(temp);
}

FemurImplantInfo FemurImplant::getImplantInfo() const
{
    return mImplantInfo;
}

void FemurImplant::FixKneeCapPath(const std::vector<PointTypeITK>& pKneeCapPath)
{
    int tSize = pKneeCapPath.size();

    bool resultOk;
    Plane bestPlane = Plane::getBestPlane(pKneeCapPath, resultOk);
    
    if (bestPlane.getIsInit() == false)
    {
        return;
    }
    
    for (int i = 0; i < tSize; i++)
    {
        Point pnt = { pKneeCapPath[i][0], pKneeCapPath[i][1], pKneeCapPath[i][2] };
        mKneeCapPath.push_back(bestPlane.getProjectionPoint(pnt));
    }
    
}

std::vector<Point> FemurImplant::GetKneeCapPath() const
{
    return mKneeCapPath;
}

Plane FemurImplant::getPlaneA() const
{
	return planeA;
}

Plane FemurImplant::getPlaneB() const
{
	return planeB;
}

Plane FemurImplant::getPlaneC() const
{
	return planeC;
}

Plane FemurImplant::getPlaneD() const
{
	return planeD;
}

Plane FemurImplant::getPlaneE() const
{
	return planeE;
}

Plane FemurImplant::getMidPlane() const
{
	return midPlane;
}

Point FemurImplant::getCortexPoint() const
{
	return cortexPoint;
}

Point FemurImplant::getPointP3() const
{
    return P3;
}

Point FemurImplant::getPointP4() const
{
    return P4;
}

Point FemurImplant::getMidPoint() const
{
	return (P3 + P4) / 2;
}

cv::Mat FemurImplant::getCortexPointMat() const
{
	cv::Mat mat(3, 1, CV_64FC1);
	mat.at <double>(0, 0) = cortexPoint.x;
	mat.at <double>(1, 0) = cortexPoint.y;
	mat.at <double>(2, 0) = cortexPoint.z;
	return mat;
}

Point FemurImplant::getDirectVectorFemurAxis() const
{
	Point directVector = planeC.getNormalVector();
	double squareNorm = directVector.x * directVector.x + directVector.y * directVector.y + directVector.z * directVector.z;
	return directVector / sqrt(squareNorm);
}

Point FemurImplant::getDirectVectorTEA() const
{
    Point directVector = (getDirectVectorFemurAxis()).cross(getDirectVectorAP());//P4 - P3;
	double squareNorm = directVector.x * directVector.x + directVector.y * directVector.y + directVector.z * directVector.z;
	return directVector / sqrt(squareNorm);
}

Point FemurImplant::getDirectVectorAP() const
{
	Point midPoint = (P3 + P4) / 2.0;
	Point tDirectVector = planeC.getNormalVector().cross(midPlane.getNormalVector());
    Plane helpPlane;
    helpPlane.init(tDirectVector, midPoint);
    Point directVector = cortexPoint - helpPlane.getProjectionPoint(cortexPoint);// Line::getFixDirectVector(tDirectVector, midPoint, cortexPoint);
	double squareNorm = directVector.x * directVector.x + directVector.y * directVector.y + directVector.z * directVector.z;
	return directVector / sqrt(squareNorm);
}

vtkSmartPointer<vtkPolyData> FemurImplant::GetImplantModel() const
{
    return implantModel;
}

vtkSmartPointer<vtkPolyData> FemurImplant::GetTransformImplantModel(const itk::Rigid3DTransform<>::Pointer transform) const
{
    itk::Matrix< double, 3, 3 > rotation = transform->GetMatrix();
    itk::Vector< double, 3 > translate = transform->GetOffset();

    vtkNew<vtkMatrix4x4> m;
    m->SetElement(0, 0, rotation[0][0]);
    m->SetElement(0, 1, rotation[0][1]);
    m->SetElement(0, 2, rotation[0][2]);

    m->SetElement(1, 0, rotation[1][0]);
    m->SetElement(1, 1, rotation[1][1]);
    m->SetElement(1, 2, rotation[1][2]);

    m->SetElement(2, 0, rotation[2][0]);
    m->SetElement(2, 1, rotation[2][1]);
    m->SetElement(2, 2, rotation[2][2]);

    m->SetElement(3, 0, 0);
    m->SetElement(3, 1, 0);
    m->SetElement(3, 2, 0);
    m->SetElement(3, 3, 1);

    m->SetElement(0, 3, translate[0]);
    m->SetElement(1, 3, translate[1]);
    m->SetElement(2, 3, translate[2]);

    vtkNew<vtkTransform> vtkTransform;
    vtkTransform->SetMatrix(m);

    vtkNew<vtkTransformPolyDataFilter> transformFilter;
    transformFilter->SetInputData(implantModel);
    transformFilter->SetTransform(vtkTransform);
    transformFilter->Update();
    return transformFilter->GetOutput();
}