
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
#include <iostream>
#include "FemurImplant.hpp"
#include "ImplantsException.hpp"


using namespace UKA::IMPLANTS;

FemurImplant::FemurImplant()
{
	isInit = false;
}

void FemurImplant::init(const Plane& pPosterior, const vtkSmartPointer<vtkPolyData> implantModel,
	const FemurImplantInfo& pImplantInfo)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_FEMUR_IMPLANT;
    }

    this->mPosterior = pPosterior;
    this->mImplantModel = implantModel;
    this->mImplantInfo = pImplantInfo;
    isInit = true;
}

FemurImplantInfo FemurImplant::getImplantInfo() const
{
    return mImplantInfo;
}

Plane FemurImplant::getPosterior() const
{
	return mPosterior;
}

vtkSmartPointer<vtkPolyData> FemurImplant::GetImplantModel() const
{
    return mImplantModel;
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
    transformFilter->SetInputData(mImplantModel);
    transformFilter->SetTransform(vtkTransform);
    transformFilter->Update();
    return transformFilter->GetOutput();
}