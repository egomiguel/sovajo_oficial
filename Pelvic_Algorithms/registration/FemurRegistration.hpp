#ifndef FEMUR_REGISTRATION_H
#define FEMUR_REGISTRATION_H

#include "Registration.hpp"
#include "registration_export.h"

class REGISTRATION_EXPORT FemurRegistration: public Registration
{
public:
    FemurRegistration(const vtkSmartPointer<vtkPolyData> img, const PointTypeITK& pHipCenterCT, const PointTypeITK& pKneeCenterCT, const PointTypeITK& pMedialEpicondyleCT);
    
    ~FemurRegistration();

    bool MakeRegistration(const std::vector<itk::Point<double, 3>>& pBonePoints, const PointTypeITK& pHipCamera, const PointTypeITK& pKneeCenterCamera, const PointTypeITK& pMedialEpicondyleCamera, bool useRandomAlignment = false);

private:
    PointTypeITK hipCenterCT;
	PointTypeITK kneeCenterCT;
	PointTypeITK medialEpicondyleCT;
};

#endif
