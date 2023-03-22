//#include <itkImageFileReader.h>
//#include <itkDirectory.h>
//#include <itkGDCMImageIO.h>
//#include "itkImageToVTKImageFilter.h"
#include "vtkFlyingEdges3D.h"
#include <itkExtractImageFilter.h>
#include <itkBinaryMask3DMeshSource.h>
#include <opencv2/opencv.hpp>
#include "Knee.hpp"
#include "ImplantsException.hpp"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkCleanPolyData.h"
#include "ImplantTools.hpp"
#include "vtkImplicitPolyDataDistance.h"
#include "TemplatePoints.hpp"
#include "LeastSquaresScaleICP.hpp"

const Point UnitVectorSagital(1, 0, 0);
const Point UnitVectorCoronal(0, 1, 0);
const Point UnitVectorAxial(0, 0, 1);

Knee::Knee()
{
    isInit = false;
}

void Knee::init(const Point& hipCenter, const Point& anteriorCortex, const Point& femurKneeCenter, const Point& lateralEpicondyle,
    const Point& medialEpicondyle, const Point& lateralCondyle, const Point& medialCondyle, const Point& lateralPlateau,
    const Point& medialPlateau, const Point& tibiaKneeCenter, const Point& tibiaTubercle, const Point& pclCenter, const Point& ankleCenter,
    const Patella& pPatella, const vtkSmartPointer<vtkPolyData> femurPoly, const vtkSmartPointer<vtkPolyData> tibiaPoly,
    double cartilage, uint8_t imageValueMax)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_KNEE;
    }

    vtkNew<vtkPolyDataConnectivityFilter> femurConnectivityFilter;
    femurConnectivityFilter->SetInputData(femurPoly);
    femurConnectivityFilter->SetExtractionModeToLargestRegion();
    femurConnectivityFilter->Update();

    vtkNew<vtkCleanPolyData> CleanFemur, CleanTibia;
    CleanFemur->SetInputData(femurConnectivityFilter->GetOutput());
    CleanFemur->Update();

    vtkNew<vtkPolyDataConnectivityFilter> tibiaConnectivityFilter;
    tibiaConnectivityFilter->SetInputData(tibiaPoly);
    tibiaConnectivityFilter->SetExtractionModeToLargestRegion();
    tibiaConnectivityFilter->Update();

    CleanTibia->SetInputData(tibiaConnectivityFilter->GetOutput());
    CleanTibia->Update();

    Plane epiHelp;
    epiHelp.init((hipCenter - femurKneeCenter), lateralEpicondyle);

    this->hipCenter = hipCenter;
    this->femurKneeCenter = femurKneeCenter;
    this->lateralEpicondyle = lateralEpicondyle;
    this->medialEpicondyle = medialEpicondyle;
    this->medialEpicondylePerp = epiHelp.getProjectionPoint(medialEpicondyle);
    this->ankleCenter = ankleCenter;
    this->lateralCondyle = lateralCondyle;
    this->medialCondyle = medialCondyle;
    this->lateralPlateau = lateralPlateau;
    this->medialPlateau = medialPlateau;
    this->pclCenter = pclCenter;
    this->femurPoly = CleanFemur->GetOutput();
    this->tibiaPoly = CleanTibia->GetOutput();
    this->femurCartilage = cartilage;
    this->tibiaCartilage = cartilage;
    this->mPatella = pPatella;
    this->anteriorCortex = anteriorCortex;
    this->tibiaKneeCenter = tibiaKneeCenter;
    this->tibiaTubercle = tibiaTubercle;

    FillFemurPoints();
    FillTibiaPoints();

    Point inferiorLateralPoint, inferiorMedialPoint;
    getDistalInferiorCondyle(inferiorLateralPoint, inferiorMedialPoint);
    this->lateralInferiorFemurPoint = inferiorLateralPoint;
    this->medialInferiorFemurPoint = inferiorMedialPoint;

    makeKneeGroovePath();

    goodSide = getGoodSide(hipCenter, femurKneeCenter, lateralEpicondyle, medialEpicondylePerp, ankleCenter);

    if (goodSide == MedialSide)
    {
        isVarus = false;
    }
    else
    {
        isVarus = true;
    }
    Point femurDirectVector = hipCenter - femurKneeCenter;
    Point directVectorEpicondyle = lateralEpicondyle - medialEpicondylePerp;
    Point tFemurDirectVectorAP = femurDirectVector.cross(directVectorEpicondyle);

    Line tempAP(tFemurDirectVectorAP, femurKneeCenter);
    Point condyleProj = tempAP.getProjectPoint(medialCondyle);

    this->femurDirectVectorAP = femurKneeCenter - condyleProj;
    findTibiaPlaneNormalVector();

    Point rightVector = medialEpicondylePerp - lateralEpicondyle;
    Point crossVector = getDirectVectorFemurAxis().cross(femurDirectVectorAP);

    double angle = ImplantTools::getAngleBetweenVectorsDegree(rightVector, crossVector);
    if (angle < 90)
    {
        isRight = true;
    }
    else
    {
        isRight = false;
    }

    getAutomaticPlateaus();/////////////////////////////////

    if (mPatella.getIsInit() == true)
    {
        mPatellaPlane = mPatella.getPatellaPlane(isRight, mPatellaDistalPosteriorPoint);
    }

    /*Plane axial, tibiaHelp, coronal;
    axial.init((hipCenter - femurKneeCenter), inferiorLateralPoint);
    tibiaHelp.init(tibiaNormalPlaneVector, medialPlateau);
    coronal.init(this->femurDirectVectorAP, lateralCondyle);

    femurLatMedDiffDistal = abs(axial.eval(medialInferiorFemurPoint));
    femurLatMedDiffPosterior = abs(coronal.eval(medialCondyle));
    tibiaLatMedDiff = abs(tibiaHelp.eval(lateralPlateau));*/

    isInit = true;
}

void Knee::init(const Point& hipCenter, const Point& anteriorCortex, const Point& femurKneeCenter, const Point& lateralEpicondyle,
    const Point& medialEpicondyle, const Point& lateralPlateau, const Point& medialPlateau, const Point& tibiaKneeCenter,
    const Point& tibiaTubercle, const Point& pclCenter, const Point& ankleCenter, const Patella& pPatella,
    const vtkSmartPointer<vtkPolyData> femurPoly, const vtkSmartPointer<vtkPolyData> tibiaPoly, KneeSideEnum pSide,
    double cartilage, uint8_t imageValueMax)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_KNEE;
    }

    vtkNew<vtkPolyDataConnectivityFilter> femurConnectivityFilter;
    femurConnectivityFilter->SetInputData(femurPoly);
    femurConnectivityFilter->SetExtractionModeToLargestRegion();
    femurConnectivityFilter->Update();

    vtkNew<vtkCleanPolyData> CleanFemur, CleanTibia;
    CleanFemur->SetInputData(femurConnectivityFilter->GetOutput());
    CleanFemur->Update();

    vtkNew<vtkPolyDataConnectivityFilter> tibiaConnectivityFilter;
    tibiaConnectivityFilter->SetInputData(tibiaPoly);
    tibiaConnectivityFilter->SetExtractionModeToLargestRegion();
    tibiaConnectivityFilter->Update();

    CleanTibia->SetInputData(tibiaConnectivityFilter->GetOutput());
    CleanTibia->Update();

    Plane epiHelp;
    epiHelp.init((hipCenter - femurKneeCenter), lateralEpicondyle);

    this->hipCenter = hipCenter;
    this->femurKneeCenter = femurKneeCenter;
    this->lateralEpicondyle = lateralEpicondyle;
    this->medialEpicondyle = medialEpicondyle;
    this->medialEpicondylePerp = epiHelp.getProjectionPoint(medialEpicondyle);
    this->ankleCenter = ankleCenter;
    this->lateralPlateau = lateralPlateau;
    this->medialPlateau = medialPlateau;
    this->pclCenter = pclCenter;
    this->femurPoly = CleanFemur->GetOutput();
    this->tibiaPoly = CleanTibia->GetOutput();
    this->femurCartilage = cartilage;
    this->tibiaCartilage = cartilage;
    this->mPatella = pPatella;
    this->anteriorCortex = anteriorCortex;
    this->tibiaKneeCenter = tibiaKneeCenter;
    this->tibiaTubercle = tibiaTubercle;

    Point forceLineVector = hipCenter - femurKneeCenter;
    Point TEA = lateralEpicondyle - medialEpicondylePerp;

    if (pSide == KneeSideEnum::KRight)
    {
        isRight = true;
        this->femurDirectVectorAP = forceLineVector.cross(TEA);
    }
    else
    {
        isRight = false;
        this->femurDirectVectorAP = TEA.cross(forceLineVector);
    }

    getAutomaticPlateaus();//////////////////////

    FillFemurPointsAndCondyles();
    FillTibiaPoints();

    makeKneeGroovePath();

    goodSide = getGoodSide(hipCenter, femurKneeCenter, lateralEpicondyle, medialEpicondylePerp, ankleCenter);

    if (goodSide == MedialSide)
    {
        isVarus = false;
    }
    else
    {
        isVarus = true;
    }

    findTibiaPlaneNormalVector();

    if (mPatella.getIsInit() == true)
    {
        mPatellaPlane = mPatella.getPatellaPlane(isRight, mPatellaDistalPosteriorPoint);
    }

    isInit = true;
}

void Knee::init(const Point& hipCenter, const Point& anteriorCortex, const Point& femurKneeCenter, const Point& lateralEpicondyle,
    const Point& medialEpicondyle, const Point& tibiaKneeCenter, const Point& tibiaTubercle, const Point& pclCenter, 
    const Point& ankleCenter, const Patella& pPatella, const vtkSmartPointer<vtkPolyData> femurPoly, 
    const vtkSmartPointer<vtkPolyData> tibiaPoly, KneeSideEnum pSide, double cartilage, uint8_t imageValueMax)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_KNEE;
    }

    vtkNew<vtkPolyDataConnectivityFilter> femurConnectivityFilter;
    femurConnectivityFilter->SetInputData(femurPoly);
    femurConnectivityFilter->SetExtractionModeToLargestRegion();
    femurConnectivityFilter->Update();

    vtkNew<vtkCleanPolyData> CleanFemur, CleanTibia;
    CleanFemur->SetInputData(femurConnectivityFilter->GetOutput());
    CleanFemur->Update();

    vtkNew<vtkPolyDataConnectivityFilter> tibiaConnectivityFilter;
    tibiaConnectivityFilter->SetInputData(tibiaPoly);
    tibiaConnectivityFilter->SetExtractionModeToLargestRegion();
    tibiaConnectivityFilter->Update();

    CleanTibia->SetInputData(tibiaConnectivityFilter->GetOutput());
    CleanTibia->Update();

    Plane epiHelp;
    epiHelp.init((hipCenter - femurKneeCenter), lateralEpicondyle);

    this->hipCenter = hipCenter;
    this->femurKneeCenter = femurKneeCenter;
    this->lateralEpicondyle = lateralEpicondyle;
    this->medialEpicondyle = medialEpicondyle;
    this->medialEpicondylePerp = epiHelp.getProjectionPoint(medialEpicondyle);
    this->ankleCenter = ankleCenter;
    this->pclCenter = pclCenter;
	this->femurPoly = CleanFemur->GetOutput();
    this->tibiaPoly = CleanTibia->GetOutput();
    this->femurCartilage = cartilage;
    this->tibiaCartilage = cartilage;
    this->mPatella = pPatella;
    this->anteriorCortex = anteriorCortex;
    this->tibiaKneeCenter = tibiaKneeCenter;
    this->tibiaTubercle = tibiaTubercle;

    Point forceLineVector = hipCenter - femurKneeCenter;
    Point TEA = lateralEpicondyle - medialEpicondylePerp;

    if (pSide == KneeSideEnum::KRight)
    {
        isRight = true;
        this->femurDirectVectorAP = forceLineVector.cross(TEA);
    }
    else
    {
        isRight = false;
        this->femurDirectVectorAP = TEA.cross(forceLineVector);
    }

    getAutomaticPlateaus();

    FillFemurPointsAndCondyles();
    FillTibiaPoints();

    //makeKneeGroovePath();

    goodSide = getGoodSide(hipCenter, femurKneeCenter, lateralEpicondyle, medialEpicondylePerp, ankleCenter);

    if (goodSide == MedialSide)
    {
        isVarus = false;
    }
    else
    {
        isVarus = true;
    }

    findTibiaPlaneNormalVector();

    if (mPatella.getIsInit() == true)
    {
        mPatellaPlane = mPatella.getPatellaPlane(isRight, mPatellaDistalPosteriorPoint);
    }

    isInit = true;
}

//void Knee::setImplantInfo(const ImplantInfo pImplant)
//{
//    implant = pImplant;
//
//    if (isVarus == false)
//    {
//        condyleResectionThickness = implant.femurPosteriorMedialThickness;// -cartilage;
//        moveCondyle = medialCondyle;
//    }
//    else
//    {
//        condyleResectionThickness = implant.femurPosteriorLateralThickness;// -cartilage;
//        moveCondyle = lateralCondyle;
//    }
//
//    Point normaliceFemurDirectVectorAP = femurDirectVectorAP;
//    normaliceFemurDirectVectorAP.normalice();
//
//    this->moveCondyle = moveCondyle + condyleResectionThickness * normaliceFemurDirectVectorAP; //getPointAtDistance(tFemurDirectVectorAP, moveCondyle, tempAP.getProjectPoint(femurKneeCenter), condyleResectionThickness);
//
//    if (isVarus == false)
//    {
//        plateauResectionThickness = implant.tibiaMedialThickness;// -cartilage;
//        movePlateau = medialPlateau;
//    }
//    else
//    {
//        plateauResectionThickness = implant.tibiaLateralThickness;// -cartilage;
//        movePlateau = lateralPlateau;
//    }
//
//    Point tibiaLineDirectVector = tibiaNormalPlaneVector;//ankleCenter - tibiaKneeCenter;
//    tibiaLineDirectVector.normalice();
//    this->movePlateau = movePlateau - plateauResectionThickness * tibiaLineDirectVector;//getPointAtDistance(tibiaLineDirectVector, movePlateau, ankleCenter, plateauResectionThickness);
//
//    if (isVarus == false)
//    {
//        inferiorMoveFemurPoint = medialInferiorFemurPoint;
//        distalCondyleResectionThickness = implant.femurDistalMedialThickness;// -cartilage;
//    }
//    else
//    {
//        inferiorMoveFemurPoint = lateralInferiorFemurPoint;
//        distalCondyleResectionThickness = implant.femurDistalLateralThickness;// -cartilage;
//    }
//
//    Point femurLineDirectVector = hipCenter - femurKneeCenter;
//    femurLineDirectVector.normalice();
//
//    this->inferiorMoveFemurPoint = inferiorMoveFemurPoint + distalCondyleResectionThickness * femurLineDirectVector;
//}

bool Knee::getIsRight() const
{
    return isRight;
}

void Knee::FillFemurPoints()
{
    vtkSmartPointer<vtkPoints> allPoints = femurPoly->GetPoints();
    vtkIdType tSize = allPoints->GetNumberOfPoints();
    for (vtkIdType i = 0; i < tSize; i++)
    {
        double pnt[3];
        allPoints->GetPoint(i, pnt);
        mFemur.push_back(Point(pnt[0], pnt[1], pnt[2]));
    }
}

void Knee::FillFemurPointsAndCondyles()
{
    Point forceLine = hipCenter - femurKneeCenter;
    Point tea = lateralEpicondyle - medialEpicondyle;

    Plane sagital, axial, coronal, transverse;
    sagital.init(tea, femurKneeCenter);
    sagital.reverseByPoint(lateralEpicondyle);

    axial.init(forceLine, femurKneeCenter);
    axial.reverseByPoint(hipCenter, false);

    coronal.init(femurDirectVectorAP, femurKneeCenter);
    coronal.reverse();

    transverse.init(forceLine, anteriorCortex);
    transverse.reverseByPoint(hipCenter, false);

    double distalLatDist = -1, distalMedDist = -1, posteriorLatDist = -1, posteriorMedDist = -1, temp;

    vtkSmartPointer<vtkPoints> allPoints = femurPoly->GetPoints();
    vtkIdType tSize = allPoints->GetNumberOfPoints();
    for (vtkIdType i = 0; i < tSize; i++)
    {
        double pnt[3];
        allPoints->GetPoint(i, pnt);

        Point myPoint = Point(pnt[0], pnt[1], pnt[2]);
        mFemur.push_back(myPoint);

        if (sagital.eval(myPoint) > 0)
        {
            temp = axial.eval(myPoint);

            if (temp > distalLatDist)
            {
                distalLatDist = temp;
                lateralInferiorFemurPoint = myPoint;
            }

            temp = coronal.eval(myPoint);

            if (temp > posteriorLatDist)
            {
                posteriorLatDist = temp;
                lateralCondyle = myPoint;
            }
        }
        else
        {
            temp = axial.eval(myPoint);

            if (temp > distalMedDist)
            {
                distalMedDist = temp;
                medialInferiorFemurPoint = myPoint;
            }

            temp = coronal.eval(myPoint);

            if (temp > posteriorMedDist)
            {
                posteriorMedDist = temp;
                medialCondyle = myPoint;
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////

    Plane coronalBase = axial.getPerpendicularPlane(lateralCondyle, medialCondyle);
    Point refPoint = medialCondyle + 1000 * femurDirectVectorAP;
    coronalBase.reverseByPoint(refPoint);

    double coronalLatDistance = -1, coronalMedDistance = -1;

    auto it1 = mFemur.begin();
    auto it2 = mFemur.end();

    for ( ; it1 != it2; ++it1 )
    {
        if (coronalBase.eval(*it1) > 0 && transverse.eval(*it1) > 0)
        {
            if (sagital.eval(*it1) > 2)
            {
                temp = coronalBase.eval(*it1);

                if (temp > coronalLatDistance)
                {
                    coronalLatDistance = temp;
                    coronalDistalLat = (*it1);
                }
            }
            else if (sagital.eval(*it1) < -2)
            {
                temp = coronalBase.eval(*it1);

                if (temp > coronalMedDistance)
                {
                    coronalMedDistance = temp;
                    coronalDistalMed = (*it1);
                }
            }

        }
    }

    //UpdateTopPointOnGroove();********************************************************************************************************
}

void Knee::UpdateTopPointOnGroove()
{
    Point forceLine = hipCenter - femurKneeCenter;
    Point tea = lateralEpicondyle - medialEpicondyle;

    Plane coronal, transverse;

    coronal.init(femurDirectVectorAP, femurKneeCenter);
    transverse.init(forceLine, anteriorCortex);

    /////////////////////////////////////////////////////////////////////////////////////////////////

    Plane lateralPlane, medialPlane, cutPlane;
    medialPlane.init(tea, coronalDistalMed);
    lateralPlane.init(tea, coronalDistalLat);
    cutPlane = coronal.getPerpendicularPlane(coronalDistalLat, coronalDistalMed);

    if (lateralPlane.eval(medialEpicondylePerp) < 0)
    {
        lateralPlane.reverse();
    }

    if (medialPlane.eval(lateralEpicondyle) < 0)
    {
        medialPlane.reverse();
    }

    vtkSmartPointer<vtkPolyData> contour = ImplantTools::getMaxContour(femurPoly, cutPlane.getNormalVector(), cutPlane.getPoint());

    ////////////////////////////////////////////////////////////////

    Point vectorAP = femurDirectVectorAP;
    vectorAP.normalice();

    vtkNew<vtkImplicitPolyDataDistance> polyDistance;
    polyDistance->SetInput(femurPoly);

    Point refPoint = (coronalDistalLat + coronalDistalMed) / 2.;
    Point extPoint = refPoint + 10. * vectorAP;
    Point midGroovePoint;

    ImplantTools::GetInterceptionWithLine(polyDistance, refPoint, extPoint, midGroovePoint);

    vtkIdType refMid = ImplantTools::GetNearestPoints(contour, midGroovePoint);
    vtkIdType refLat = ImplantTools::GetNearestPoints(contour, coronalDistalLat);
    vtkIdType refMed = ImplantTools::GetNearestPoints(contour, coronalDistalMed);

    std::list<std::pair<vtkIdType, vtkIdType>> lines;
    ImplantTools::ExtractSortLines(contour, lines);

    int cont = lines.size();

    while (cont > -1 && lines.front().first != refMid)
    {
        lines.push_back(lines.front());
        lines.pop_front();
        cont--;
    }

    if (cont < 0)
    {
        throw ImplantExceptionCode::CHECK_LANDMARKS_CAN_NOT_DETERMINE_BEGIN_POINT_OF_KNEE_GROOVE;
    }

    std::vector<Point> curve;

    auto it1 = lines.begin();
    auto it2 = lines.end();

    for (; it1 != it2; ++it1)
    {
        if (it1->first != refLat && it1->first != refMed)
        {
            double pnt[3];
            contour->GetPoint((it1)->first, pnt);
            curve.push_back(Point(pnt[0], pnt[1], pnt[2]));
        }
        else
        {
            break;
        }
    }

    auto rit1 = lines.rbegin();
    auto rit2 = lines.rend();

    for (; rit1 != rit2; ++rit1)
    {
        if (rit1->first != refLat && rit1->first != refMed)
        {
            double pnt[3];
            contour->GetPoint((rit1)->first, pnt);
            curve.push_back(Point(pnt[0], pnt[1], pnt[2]));
        }
        else
        {
            break;
        }
    }

    ///////////////////////////////////////////////////////////////

    Plane basePlane = transverse.getPerpendicularPlane(lateralCondyle, medialCondyle);
    basePlane.reverseByPoint(femurKneeCenter);

    Point topPointTemp;

    try
    {
        topPointTemp = ImplantTools::getLocalMinimum(curve, cutPlane, cutPlane.getProjectionVector(basePlane.getNormalVector()));
    }
    catch (...)
    {
        throw ImplantExceptionCode::CHECK_LANDMARKS_CAN_NOT_DETERMINE_BEGIN_POINT_OF_KNEE_GROOVE;
    }

    double myClosest[3];
    double pnt[3] = { topPointTemp.x, topPointTemp.y, topPointTemp.z };

    polyDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);
    topPointOnPatellaPath = Point(myClosest[0], myClosest[1], myClosest[2]);
}

void Knee::FillTibiaPoints()
{
    vtkSmartPointer<vtkPoints> allPoints = tibiaPoly->GetPoints();
    vtkIdType tSize = allPoints->GetNumberOfPoints();
    for (vtkIdType i = 0; i < tSize; i++)
    {
        double pnt[3];
        allPoints->GetPoint(i, pnt);
        mTibia.push_back(Point(pnt[0], pnt[1], pnt[2]));
    }
}

Point Knee::getPatellaRotationPoint() const
{
    return mPatella.getPatellaRotationPoint();
}

void Knee::setTibiaSlope(double angleInDegrees)
{
    findTibiaPlaneNormalVector(angleInDegrees);
}

void Knee::findTibiaPlaneNormalVector(double degreeAngle)
{
    Point direct_vector_axis = ankleCenter - tibiaKneeCenter;
    direct_vector_axis.normalice();
    Plane tibia_help;
    tibia_help.init(direct_vector_axis, tibiaKneeCenter);
    Point pseudo_vector_AP;
    if (degreeAngle >= 0)
    {
        pseudo_vector_AP = pclCenter - tibiaTubercle;
    }
    else
    {
        pseudo_vector_AP = tibiaTubercle - pclCenter;
        degreeAngle = 0 - degreeAngle;
    }
    Point vector_AP_respect_axis = tibia_help.getProjectionVector(pseudo_vector_AP);
    vector_AP_respect_axis.normalice();

    double tibiaSize = Line::getDistanceBetweenPoints(tibiaKneeCenter, ankleCenter);
    double complementAngleRadians = ((90.0 - degreeAngle) * PI) / 180.0;
    double angleRadians = (degreeAngle * PI) / 180.0;
    double lenghtRotateVector = tibiaSize * cos(complementAngleRadians);
    double lenghtSideAP = lenghtRotateVector * cos(angleRadians);
    double lenghtParallelAxis = lenghtRotateVector * sin(angleRadians);

    Point point_AP = tibiaKneeCenter + lenghtSideAP * vector_AP_respect_axis;
    tibiaRotatePoint = point_AP + lenghtParallelAxis * direct_vector_axis;
    tibiaNormalPlaneVector = tibiaRotatePoint - ankleCenter;
}

/*
void Knee::findTibiaDirectVectorAP()
{
    Point directVector = tibiaTubercle - pclCenter;
    Point tibiaAxisVector = (tibiaKneeCenter - ankleCenter);
    Plane tibiaHelpPlane;
    tibiaHelpPlane.init(tibiaAxisVector, tibiaKneeCenter);
    Point directVectorAP = tibiaHelpPlane.getProjectionVector(directVector);

    tibiaFixedDirectVectorAP = directVectorAP;//Line::getFixDirectVector(directVectorAP, tibiaKneeCenter, tibiaTubercle);
}*/
/*
void Knee::readImagePoints(const ImplantImageType::Pointer& image, uint8_t imageValueMax, std::vector<Point>& meshPoints) const
{
    using MeshType = itk::Mesh<double, 3>;
    auto meshSource = itk::BinaryMask3DMeshSource<ImplantImageType, MeshType>::New();
    meshSource->SetObjectValue(imageValueMax);
    meshSource->SetInput(image);

    try
    {
        meshSource->Update();
    }
    catch (itk::ExceptionObject & exp)
    {
        std::string msg = exp.what();
        throw ImplantsException(msg);
    }

    itk::Mesh<double, 3>::ConstPointer mesh = meshSource->GetOutput();
    MeshType::PointsContainer::ConstPointer points;
    points = mesh->GetPoints();

    if (points->size() == 0)
    {
        throw ImplantsException("Could not get the list of surface points.");
    }

    itk::Point<double, 3> myPoint;
    auto pointsBegin = points->Begin();
    auto pointsEnd = points->End();
    while (pointsBegin != pointsEnd)
    {
        myPoint = pointsBegin.Value();
        meshPoints.push_back(Point(myPoint[0], myPoint[1], myPoint[2]));
        ++pointsBegin;
    }
}
*/
Point Knee::getFemurVectorTEA() const
{
    Point directVector = (getDirectVectorFemurAxis()).cross(getFemurDirectVectorAP());
    double squareNorm = directVector.x * directVector.x + directVector.y * directVector.y + directVector.z * directVector.z;
    return directVector / sqrt(squareNorm);
}

Point Knee::getTibiaVectorTEA() const
{
    Point directVector = tibiaNormalPlaneVector.cross(getTibiaDirectVectorAP());
    double squareNorm = directVector.x * directVector.x + directVector.y * directVector.y + directVector.z * directVector.z;
    return directVector / sqrt(squareNorm);
}

Point Knee::getFemurVectorLateralTEA() const
{
    Point vector = lateralEpicondyle - medialEpicondylePerp;
    vector.normalice();
    return vector;
}

Point Knee::getTibiaVectorLateralTEA() const
{
    Point vector = tibiaNormalPlaneVector.cross(getTibiaDirectVectorAP());
    vector.normalice();
    if (isRight == true)
    {
        return -vector;
    }
    else
    {
        return vector;
    }

}

Point Knee::getDirectVectorFemurAxis() const
{
    Point directVector = hipCenter - femurKneeCenter;
    double squareNorm = directVector.x * directVector.x + directVector.y * directVector.y + directVector.z * directVector.z;
    return directVector / sqrt(squareNorm);
}

Point Knee::getFemurDirectVectorAP() const
{
    Point directVector = femurDirectVectorAP;
    double squareNorm = directVector.x * directVector.x + directVector.y * directVector.y + directVector.z * directVector.z;
    return directVector / sqrt(squareNorm);
}

Point Knee::getTibiaDirectVectorAP() const
{
    Plane tibiaPlane;
    tibiaPlane.init(tibiaNormalPlaneVector, medialPlateau);

    Point directVectorTemp = tibiaTubercle - pclCenter;
    Point directVector = tibiaPlane.getProjectionVector(directVectorTemp);

    double squareNorm = directVector.x * directVector.x + directVector.y * directVector.y + directVector.z * directVector.z;
    return directVector / sqrt(squareNorm);
}

Point Knee::getNormalVectorTibiaPlane() const
{
    Point directVector = tibiaNormalPlaneVector;
    double squareNorm = directVector.x * directVector.x + directVector.y * directVector.y + directVector.z * directVector.z;
    return directVector / sqrt(squareNorm);
}

Point Knee::getAnteriorCortex() const
{
    return anteriorCortex;
}

//Point Knee::getMoveCondyle1() const
//{
//    return moveCondyle;
//}

Point Knee::getMoveCondyle(const FemurImplantInfo& pImplant) const
{
    double resectionThickness;
    Point myMoveCondyle;

    /*
    if (goodSide == MedialSide)
    {
        //resectionThickness = pImplant.femurPosteriorMedialThickness + LAT_MED_DIFFERENCE - cartilage;
        resectionThickness = pImplant.femurPosteriorMedialThickness - femurCartilage;
        myMoveCondyle = medialCondyle;
    }
    else
    {
        //resectionThickness = pImplant.femurPosteriorLateralThickness - cartilage;
        resectionThickness = pImplant.femurPosteriorLateralThickness - femurLatMedDiffPosterior - femurCartilage;

        //if (resectionThickness < 2)
        //{
        //    resectionThickness = 2;
        //}

        myMoveCondyle = lateralCondyle;
    }
    */

    Point normaliceFemurDirectVectorAP = femurDirectVectorAP;
    normaliceFemurDirectVectorAP.normalice();

    ///////////////////////////////////////////////////////////////////
    Plane helpPlane;
    helpPlane.init(normaliceFemurDirectVectorAP, femurKneeCenter);

    if (abs(helpPlane.eval(medialCondyle)) > abs(helpPlane.eval(lateralCondyle)))
    {
        myMoveCondyle = medialCondyle;
        resectionThickness = pImplant.femurPosteriorMedialThickness - femurCartilage;
    }
    else
    {
        myMoveCondyle = lateralCondyle;
        resectionThickness = pImplant.femurPosteriorLateralThickness - femurCartilage;
    }
    //////////////////////////////////////////////////////////////////////

    myMoveCondyle = myMoveCondyle + resectionThickness * normaliceFemurDirectVectorAP;

    return myMoveCondyle;
}

Point Knee::getFemurKneeCenter() const
{
    return femurKneeCenter;
}

Point Knee::getTibiaKneeCenter() const
{
    return tibiaKneeCenter;
}

Point Knee::getLateralCondyle() const
{
    return lateralCondyle;
}

Point Knee::getMedialCondyle() const
{
    return medialCondyle;
}

Point Knee::getLateralPlateau() const
{
    return lateralPlateau;
}

Point Knee::getMedialPlateau() const
{
    return medialPlateau;
}

//Point Knee::getMovePlateau1() const
//{
//    return movePlateau;
//}

Point Knee::getMovePlateau(const TibiaImplantInfo& pImplant) const
{
    double resectionThickness;
    Point myMovePlateau;
    /*
    if (goodSide == MedialSide)
    {
        //resectionThickness = pImplant.tibiaMedialThickness - LAT_MED_DIFFERENCE - cartilage;
        resectionThickness = pImplant.tibiaMedialThickness - tibiaLatMedDiff - tibiaCartilage;

        //if (resectionThickness < 2)
        //{
        //    resectionThickness = 2;
        //}

        myMovePlateau = medialPlateau;
    }
    else
    {
        //resectionThickness = pImplant.tibiaLateralThickness - cartilage;
        resectionThickness = pImplant.tibiaLateralThickness - tibiaCartilage;
        myMovePlateau = lateralPlateau;
    }
    */

    Point tibiaLineDirectVector = tibiaNormalPlaneVector;
    tibiaLineDirectVector.normalice();

    ///////////////////////////////////////////////////////////////////
    Plane helpPlane;
    helpPlane.init(tibiaLineDirectVector, tibiaTubercle);

    if (abs(helpPlane.eval(medialPlateau)) > abs(helpPlane.eval(lateralPlateau)))
    {
        myMovePlateau = medialPlateau;
        resectionThickness = pImplant.tibiaMedialThickness - tibiaCartilage;
    }
    else
    {
        myMovePlateau = lateralPlateau;
        resectionThickness = pImplant.tibiaLateralThickness - tibiaCartilage;
    }
    //////////////////////////////////////////////////////////////////////

    myMovePlateau = myMovePlateau - resectionThickness * tibiaLineDirectVector;

    return myMovePlateau;
}

Point Knee::getAnkleCenter() const
{
    return ankleCenter;
}

Point Knee::getLateralEpicondyle() const
{
    return lateralEpicondyle;
}

Point Knee::getMedialEpicondyle() const
{
    return medialEpicondyle;
}

Point Knee::getMedialEpicondylePerp() const
{
    return medialEpicondylePerp;
}

Point Knee::getHipCenter() const
{
    return hipCenter;
}

Point Knee::getTibiaTubercle() const
{
    return tibiaTubercle;
}

Point Knee::getLateralInferiorFemurPoint() const
{
    return lateralInferiorFemurPoint;
}

Point Knee::getMedialInferiorFemurPoint() const
{
    return medialInferiorFemurPoint;
}

void Knee::setLateralAndMedialInferiorFemurPoints(const Point& pLateral, const Point& pMedial)
{
    lateralInferiorFemurPoint = pLateral;
    medialInferiorFemurPoint = pMedial;
}

void Knee::setLateralAndMedialPlateauPoints(const Point& pLateral, const Point& pMedial)
{
    this->lateralPlateau = pLateral;
    this->medialPlateau = pMedial;
}

void Knee::setLateralAndMedialPosteriorFemurPoints(const Point& pLateral, const Point& pMedial)
{
    lateralCondyle = pLateral;
    medialCondyle = pMedial;
    
    /////////////////////////////////////////////////////////////////

    Point forceLine = hipCenter - femurKneeCenter;
    Point tea = lateralEpicondyle - medialEpicondyle;

    Plane sagital, axial, transverse;

    sagital.init(tea, femurKneeCenter);
    sagital.reverseByPoint(lateralEpicondyle);

    axial.init(forceLine, femurKneeCenter);

    transverse.init(forceLine, anteriorCortex);
    transverse.reverseByPoint(hipCenter, false);

    double temp;

    Plane coronalBase = axial.getPerpendicularPlane(lateralCondyle, medialCondyle);
    Point refPoint = medialCondyle + 1000 * femurDirectVectorAP;
    coronalBase.reverseByPoint(refPoint);

    double coronalLatDistance = -1, coronalMedDistance = -1;

    auto it1 = mFemur.begin();
    auto it2 = mFemur.end();

    for (; it1 != it2; ++it1)
    {
        if (coronalBase.eval(*it1) > 0 && transverse.eval(*it1) > 0)
        {
            if (sagital.eval(*it1) > 2)
            {
                temp = coronalBase.eval(*it1);

                if (temp > coronalLatDistance)
                {
                    coronalLatDistance = temp;
                    coronalDistalLat = (*it1);
                }
            }
            else if (sagital.eval(*it1) < -2)
            {
                temp = coronalBase.eval(*it1);

                if (temp > coronalMedDistance)
                {
                    coronalMedDistance = temp;
                    coronalDistalMed = (*it1);
                }
            }

        }
    }

    //UpdateTopPointOnGroove();

    //makeKneeGroovePath();
}

//Point Knee::getInferiorMoveFemurPoint1() const
//{
//    return inferiorMoveFemurPoint;
//}

Point Knee::getInferiorMoveFemurPoint(const FemurImplantInfo& pImplant) const
{
    double resectionThickness;
    Point myInferiorMoveFemurPoint;

    /*
    if (goodSide == MedialSide)
    {
        myInferiorMoveFemurPoint = medialInferiorFemurPoint;
        //resectionThickness = pImplant.femurDistalMedialThickness + LAT_MED_DIFFERENCE - cartilage;
        resectionThickness = pImplant.femurDistalMedialThickness - femurCartilage;
    }
    else
    {
        myInferiorMoveFemurPoint = lateralInferiorFemurPoint;
        //resectionThickness = pImplant.femurDistalLateralThickness - cartilage;
        resectionThickness = pImplant.femurDistalLateralThickness - femurLatMedDiffDistal - femurCartilage;

        //if (resectionThickness < 2)
        //{
        //    resectionThickness = 2;
        //}
    }
    */

    Point femurLineDirectVector = hipCenter - femurKneeCenter;
    femurLineDirectVector.normalice();

    ///////////////////////////////////////////////////////////////////
    Plane helpPlane;
    helpPlane.init(femurLineDirectVector, anteriorCortex);

    if (abs(helpPlane.eval(medialInferiorFemurPoint)) > abs(helpPlane.eval(lateralInferiorFemurPoint)))
    {
        myInferiorMoveFemurPoint = medialInferiorFemurPoint;
        resectionThickness = pImplant.femurDistalMedialThickness - femurCartilage;
    }
    else
    {
        myInferiorMoveFemurPoint = lateralInferiorFemurPoint;
        resectionThickness = pImplant.femurDistalLateralThickness - femurCartilage;
    }
    //////////////////////////////////////////////////////////////////////

    myInferiorMoveFemurPoint = myInferiorMoveFemurPoint + resectionThickness * femurLineDirectVector;

    return myInferiorMoveFemurPoint;
}

Point Knee::getTibiaRotatePoint() const
{
    return tibiaRotatePoint;
}

Point Knee::getPclCenterPoint() const
{
    return pclCenter;
}

//Point Knee::getMoveTibiaKneeCenter1() const
//{
//    Plane tibiaHelp;
//    Point normalVector = tibiaNormalPlaneVector;//tibiaKneeCenter - ankleCenter;
//    tibiaHelp.init(normalVector, movePlateau);
//    Point moveKnee = tibiaHelp.getProjectionPoint(tibiaKneeCenter);
//    return moveKnee;
//}

Point Knee::getMoveTibiaKneeCenter(const TibiaImplantInfo& pImplant) const
{
    double resectionThickness;
    Point myMovePlateau;

    /*
    if (goodSide == MedialSide)
    {
        //resectionThickness = pImplant.tibiaMedialThickness - LAT_MED_DIFFERENCE - cartilage;
        resectionThickness = pImplant.tibiaMedialThickness - tibiaLatMedDiff - tibiaCartilage;

        //if (resectionThickness < 2)
        //{
         //   resectionThickness = 2;
        //}

        myMovePlateau = medialPlateau;
    }
    else
    {
        //resectionThickness = pImplant.tibiaLateralThickness - cartilage;
        resectionThickness = pImplant.tibiaLateralThickness - tibiaCartilage;
        myMovePlateau = lateralPlateau;
    }
    */

    Point tibiaLineDirectVector = tibiaNormalPlaneVector;
    tibiaLineDirectVector.normalice();

    ///////////////////////////////////////////////////////////////////
    Plane helpPlane;
    helpPlane.init(tibiaLineDirectVector, tibiaTubercle);

    if (abs(helpPlane.eval(medialPlateau)) > abs(helpPlane.eval(lateralPlateau)))
    {
        myMovePlateau = medialPlateau;
        resectionThickness = pImplant.tibiaMedialThickness - tibiaCartilage;
    }
    else
    {
        myMovePlateau = lateralPlateau;
        resectionThickness = pImplant.tibiaLateralThickness - tibiaCartilage;
    }
    //////////////////////////////////////////////////////////////////////


    myMovePlateau = myMovePlateau - resectionThickness * tibiaLineDirectVector;

    Plane tibiaHelp;
    tibiaHelp.init(tibiaLineDirectVector, myMovePlateau);
    Point moveKnee = tibiaHelp.getProjectionPoint(tibiaKneeCenter);
    return moveKnee;
}

cv::Mat Knee::getTibiaCenterPointOnImplantAP(const TibiaImplantInfo& pImplant) const
{
    Point moveKnee = getMoveTibiaKneeCenter(pImplant);
    Plane tibiaCutPlane;
    tibiaCutPlane.init(tibiaNormalPlaneVector, moveKnee);
    //tibiaCutPlane.normalizeNormalVector();
    Point tibiaTubercleProj = tibiaCutPlane.getProjectionPoint(tibiaTubercle);
    Point pclCenterProj = tibiaCutPlane.getProjectionPoint(pclCenter);

    Point tibiaVectorAP = tibiaTubercleProj - pclCenterProj;
    Line tibiaAP(tibiaVectorAP, tibiaTubercleProj);
    Point normalXY(0.0, 0.0, 1.0);
    Point rotationAxis = normalXY.cross(tibiaCutPlane.getNormalVector());
    rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
    double rotationAngle = Line::getAngleBetweenVectors(normalXY, tibiaCutPlane.getNormalVector());
    cv::Mat rotation_1 = Line::getRotateMatrix(rotationAxis, -rotationAngle);
    cv::Mat rotation_2 = Line::getRotateMatrix(rotationAxis, rotationAngle);
    cv::Mat rotateVector_1 = rotation_1 * tibiaCutPlane.getNormalVectorMat();
    cv::Mat rotateVector_2 = rotation_2 * tibiaCutPlane.getNormalVectorMat();

    cv::Mat rotate;

    double distance_1 = Line::getAngleBetweenVectors(Point(rotateVector_1), normalXY);
    double distance_2 = Line::getAngleBetweenVectors(Point(rotateVector_2), normalXY);

    if (distance_1 < distance_2)
    {
        rotate = rotation_1;
    }
    else
    {
        rotate = rotation_2;
    }
    std::vector<cv::Point2f> coplanar2d;

    //////////////////////////////////////////////

    vtkSmartPointer<vtkPolyData> contour = ImplantTools::getMaxContour(tibiaPoly, tibiaCutPlane.getNormalVector(), tibiaCutPlane.getPoint());
    vtkSmartPointer<vtkPoints> contourPoints = contour->GetPoints();
    vtkIdType tSize = contourPoints->GetNumberOfPoints();

    double x, y;
    double meanZ = 0.0;
    for (vtkIdType i = 0; i < tSize; i++)
    {
        double pnt[3];
        contourPoints->GetPoint(i, pnt);

        Point pointProj(pnt[0], pnt[1], pnt[2]);
        cv::Mat point(3, 1, CV_64F);
        point.at <double>(0, 0) = pointProj.x;
        point.at <double>(1, 0) = pointProj.y;
        point.at <double>(2, 0) = pointProj.z;
        cv::Mat rotatePointMat = rotate * point;
        x = rotatePointMat.at <double>(0, 0);
        y = rotatePointMat.at <double>(1, 0);
        coplanar2d.push_back(cv::Point2f(float(x), float(y)));
        meanZ = meanZ + rotatePointMat.at <double>(2, 0);
    }

    //////////////////////////////////////////////

    /*auto it1 = mTibia.begin();
    auto it2 = mTibia.end();
    double x, y;
    double meanZ = 0.0;
    for (; it1 != it2; ++it1)
    {
        if (tibiaCutPlane.isPointNearToPlane(*it1, 1.0) == true)
        {
            Point pointProj = tibiaCutPlane.getProjectionPoint(*it1);
            cv::Mat point(3, 1, CV_64F);
            point.at <double>(0, 0) = pointProj.x;
            point.at <double>(1, 0) = pointProj.y;
            point.at <double>(2, 0) = pointProj.z;
            cv::Mat rotatePointMat = rotate * point;
            x = rotatePointMat.at <double>(0, 0);
            y = rotatePointMat.at <double>(1, 0);
            coplanar2d.push_back(cv::Point2f(float(x), float(y)));
            meanZ = meanZ + rotatePointMat.at <double>(2, 0);
        }
    }*/
    cv::Mat result(3, 1, CV_64F);
    unsigned int coplanarSize = coplanar2d.size();
    if (coplanarSize == 0)
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_MIDPOINT_OF_AP_AXIS_ON_TIBIA_IMPLANT_MATCH;
        //result.at <double>(0, 0) = 0;
        //result.at <double>(1, 0) = 0;
        //result.at <double>(2, 0) = 0;
        //return result;
    }
    meanZ = meanZ / double(coplanarSize);

    cv::RotatedRect box = cv::fitEllipse(coplanar2d);
    float centerX, centerY;
    centerX = box.center.x;
    centerY = box.center.y;
    cv::Mat centerTemp3d(3, 1, CV_64F);
    centerTemp3d.at <double>(0, 0) = double(centerX);
    centerTemp3d.at <double>(1, 0) = double(centerY);
    centerTemp3d.at <double>(2, 0) = meanZ;

    cv::Mat centerMat3d = rotate.inv() * centerTemp3d;

    Point center3d = Point(centerMat3d);

    Point tibiaCenter = tibiaAP.getProjectPoint(center3d);
    result.at <double>(0, 0) = tibiaCenter.x;
    result.at <double>(1, 0) = tibiaCenter.y;
    result.at <double>(2, 0) = tibiaCenter.z;

    return result;
}

Plane Knee::getEquisPlaneTibia() const
{
    Point vectorAxis = tibiaNormalPlaneVector;
    vectorAxis.normalice();
    Plane tibiaCutPlane;
    double averageDist = 10.0;
    Point pointBigDistance = getLateralPlateau() - averageDist * vectorAxis;

    tibiaCutPlane.init(vectorAxis, pointBigDistance);

    Point pointSmallDistance;
    double bigDistance, smallDistance;
    bigDistance = averageDist;

    if (abs(tibiaCutPlane.eval(getMedialPlateau())) > averageDist)
    {
        pointBigDistance = getMedialPlateau() - averageDist * vectorAxis;
        tibiaCutPlane.movePlane(pointBigDistance);
        pointSmallDistance = tibiaCutPlane.getProjectionPoint(getLateralPlateau());
        smallDistance = abs(tibiaCutPlane.eval(getLateralPlateau()));
    }
    else
    {
        pointSmallDistance = tibiaCutPlane.getProjectionPoint(getMedialPlateau());
        smallDistance = abs(tibiaCutPlane.eval(getMedialPlateau()));
    }

    double moveDistUp = (bigDistance - smallDistance) / 2.0;

    if (moveDistUp < 0.05)
    {
        return tibiaCutPlane;
    }

    Point centerOnPlane = (pointBigDistance + pointSmallDistance) / 2.0;
    Point movePointUp = pointBigDistance + vectorAxis * moveDistUp;

    double bigDistanceToCenter = ImplantTools::getDistanceBetweenPoints(pointBigDistance, centerOnPlane);
    double angleUpRad = atan((moveDistUp / bigDistanceToCenter));
    double angleDownRand = (90.0 - (angleUpRad * 180.0 / PI)) * PI / 180.0;
    double moveDistDown = bigDistanceToCenter * tan(angleDownRand);
    Point newPointDown = pointBigDistance - vectorAxis * moveDistDown;
    vectorAxis = centerOnPlane - newPointDown;

    tibiaCutPlane.deletePlane();
    tibiaCutPlane.init(vectorAxis, movePointUp);

    return tibiaCutPlane;
}

std::pair<double, double> Knee::getTibiaAutomaticAxis(Point& pVectorApFront, Point& pVectorTeaLat, Point& pCenterPoint) const
{
    Plane tibiaCutPlane = getEquisPlaneTibia();

    Point normalXY(0.0, 0.0, 1.0);
    Point rotationAxis = normalXY.cross(tibiaCutPlane.getNormalVector());
    rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
    double rotationAngle = Line::getAngleBetweenVectors(normalXY, tibiaCutPlane.getNormalVector());
    cv::Mat rotation_1 = Line::getRotateMatrix(rotationAxis, -rotationAngle);
    cv::Mat rotation_2 = Line::getRotateMatrix(rotationAxis, rotationAngle);
    cv::Mat rotateVector_1 = rotation_1 * tibiaCutPlane.getNormalVectorMat();
    cv::Mat rotateVector_2 = rotation_2 * tibiaCutPlane.getNormalVectorMat();

    cv::Mat rotate;

    double distance_1 = Line::getAngleBetweenVectors(Point(rotateVector_1), normalXY);
    double distance_2 = Line::getAngleBetweenVectors(Point(rotateVector_2), normalXY);

    if (distance_1 < distance_2)
    {
        rotate = rotation_1;
    }
    else
    {
        rotate = rotation_2;
    }
    std::vector<cv::Point2f> coplanar2d;

    //////////////////////////////////////////////

    vtkSmartPointer<vtkPolyData> contour = ImplantTools::getMaxContour(tibiaPoly, tibiaCutPlane.getNormalVector(), tibiaCutPlane.getPoint());
    vtkSmartPointer<vtkPoints> contourPoints = contour->GetPoints();
    vtkIdType tSize = contourPoints->GetNumberOfPoints();

    double x, y;
    double meanZ = 0.0;
    for (vtkIdType i = 0; i < tSize; i++)
    {
        double pnt[3];
        contourPoints->GetPoint(i, pnt);

        Point pointProj(pnt[0], pnt[1], pnt[2]);
        cv::Mat point(3, 1, CV_64F);
        point.at <double>(0, 0) = pointProj.x;
        point.at <double>(1, 0) = pointProj.y;
        point.at <double>(2, 0) = pointProj.z;
        cv::Mat rotatePointMat = rotate * point;
        x = rotatePointMat.at <double>(0, 0);
        y = rotatePointMat.at <double>(1, 0);
        coplanar2d.push_back(cv::Point2f(float(x), float(y)));
        meanZ = meanZ + rotatePointMat.at <double>(2, 0);
    }

    unsigned int coplanarSize = coplanar2d.size();
    if (coplanarSize == 0)
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_ANATOMIC_AXIS_ON_TIBIA;
    }
    meanZ = meanZ / double(coplanarSize);

    cv::RotatedRect box = cv::fitEllipse(coplanar2d);
    cv::Point2f vtx[4];
    box.points(vtx);
    Point bottomLeft, topLeft, topRight, bottomRight, centerPoint;

    bottomLeft = Point(vtx[0].x, vtx[0].y, meanZ);
    topLeft = Point(vtx[1].x, vtx[1].y, meanZ);
    topRight = Point(vtx[2].x, vtx[2].y, meanZ);
    bottomRight = Point(vtx[3].x, vtx[3].y, meanZ);

    cv::Point2f ellipseCenter = box.center;
    centerPoint = Point(ellipseCenter.x, ellipseCenter.y, meanZ);

    Point AP1, AP2, TEA1, TEA2;

    if (ImplantTools::getDistanceBetweenPoints(bottomLeft, topLeft) < ImplantTools::getDistanceBetweenPoints(topLeft, topRight))
    {
        TEA1 = (bottomLeft + topLeft) / 2.0;
        TEA2 = (topRight + bottomRight) / 2.0;

        AP1 = (topRight + topLeft) / 2.0;
        AP2 = (bottomLeft + bottomRight) / 2.0;
    }
    else
    {
        TEA1 = (topRight + topLeft) / 2.0;
        TEA2 = (bottomLeft + bottomRight) / 2.0;

        AP1 = (bottomLeft + topLeft) / 2.0;
        AP2 = (topRight + bottomRight) / 2.0;
    }

    cv::Mat aMat = rotate.inv() * AP1.ToMatPoint();
    cv::Mat bMat = rotate.inv() * AP2.ToMatPoint();

    cv::Mat cMat = rotate.inv() * TEA1.ToMatPoint();
    cv::Mat dMat = rotate.inv() * TEA2.ToMatPoint();

    cv::Mat centerMat = rotate.inv() * centerPoint.ToMatPoint();

    AP1 = Point(aMat);
    AP2 = Point(bMat);
    TEA1 = Point(cMat);
    TEA2 = Point(dMat);

    Point vectorAP = AP1 - AP2;
    Point vectorTEA = TEA1 - TEA2;

    Plane coronal, sagital;
    coronal.init(vectorAP, pclCenter);
    sagital.init(vectorTEA, medialPlateau);

    coronal.reverseByPoint(tibiaTubercle);
    sagital.reverseByPoint(lateralPlateau);

    pVectorApFront = coronal.getNormalVector();
    pVectorTeaLat = sagital.getNormalVector();
    pCenterPoint = Point(centerMat);

    cv::Size2f ellipseSize = box.size;

    if (ellipseSize.width < ellipseSize.height)
    {
        return std::make_pair(ellipseSize.width, ellipseSize.height);
    }
    else
    {
        return std::make_pair(ellipseSize.height, ellipseSize.width);
    }
}

cv::Mat Knee::getTibiaCenterPointApAutomatic() const
{
    Plane tibiaCutPlane = getEquisPlaneTibia();

    Point tibiaTubercleProj = tibiaCutPlane.getProjectionPoint(tibiaTubercle);
    Point pclCenterProj = tibiaCutPlane.getProjectionPoint(pclCenter);

    Point tibiaVectorAP = tibiaTubercleProj - pclCenterProj;
    Line tibiaAP(tibiaVectorAP, tibiaTubercleProj);
    Point normalXY(0.0, 0.0, 1.0);
    Point rotationAxis = normalXY.cross(tibiaCutPlane.getNormalVector());
    rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
    double rotationAngle = Line::getAngleBetweenVectors(normalXY, tibiaCutPlane.getNormalVector());
    cv::Mat rotation_1 = Line::getRotateMatrix(rotationAxis, -rotationAngle);
    cv::Mat rotation_2 = Line::getRotateMatrix(rotationAxis, rotationAngle);
    cv::Mat rotateVector_1 = rotation_1 * tibiaCutPlane.getNormalVectorMat();
    cv::Mat rotateVector_2 = rotation_2 * tibiaCutPlane.getNormalVectorMat();

    cv::Mat rotate;

    double distance_1 = Line::getAngleBetweenVectors(Point(rotateVector_1), normalXY);
    double distance_2 = Line::getAngleBetweenVectors(Point(rotateVector_2), normalXY);

    if (distance_1 < distance_2)
    {
        rotate = rotation_1;
    }
    else
    {
        rotate = rotation_2;
    }
    std::vector<cv::Point2f> coplanar2d;

    //////////////////////////////////////////////

    vtkSmartPointer<vtkPolyData> contour = ImplantTools::getMaxContour(tibiaPoly, tibiaCutPlane.getNormalVector(), tibiaCutPlane.getPoint());
    vtkSmartPointer<vtkPoints> contourPoints = contour->GetPoints();
    vtkIdType tSize = contourPoints->GetNumberOfPoints();

    double x, y;
    double meanZ = 0.0;
    for (vtkIdType i = 0; i < tSize; i++)
    {
        double pnt[3];
        contourPoints->GetPoint(i, pnt);

        Point pointProj(pnt[0], pnt[1], pnt[2]);
        cv::Mat point(3, 1, CV_64F);
        point.at <double>(0, 0) = pointProj.x;
        point.at <double>(1, 0) = pointProj.y;
        point.at <double>(2, 0) = pointProj.z;
        cv::Mat rotatePointMat = rotate * point;
        x = rotatePointMat.at <double>(0, 0);
        y = rotatePointMat.at <double>(1, 0);
        coplanar2d.push_back(cv::Point2f(float(x), float(y)));
        meanZ = meanZ + rotatePointMat.at <double>(2, 0);
    }

    cv::Mat result(3, 1, CV_64F);
    unsigned int coplanarSize = coplanar2d.size();
    if (coplanarSize == 0)
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_MIDPOINT_OF_AP_AXIS_ON_TIBIA;
    }
    meanZ = meanZ / double(coplanarSize);

    cv::RotatedRect box = cv::fitEllipse(coplanar2d);
    float centerX, centerY;
    centerX = box.center.x;
    centerY = box.center.y;
    cv::Mat centerTemp3d(3, 1, CV_64F);
    centerTemp3d.at <double>(0, 0) = double(centerX);
    centerTemp3d.at <double>(1, 0) = double(centerY);
    centerTemp3d.at <double>(2, 0) = meanZ;

    cv::Mat centerMat3d = rotate.inv() * centerTemp3d;

    Point center3d = Point(centerMat3d);

    Point tibiaCenter = tibiaAP.getProjectPoint(center3d);
    result.at <double>(0, 0) = tibiaCenter.x;
    result.at <double>(1, 0) = tibiaCenter.y;
    result.at <double>(2, 0) = tibiaCenter.z;
    return result;
}

//Point Knee::getMoveFemurKneeCenter1() const
//{
//    Plane femurHelp;
//    Point normalVector = hipCenter - femurKneeCenter;
//    femurHelp.init(normalVector, inferiorMoveFemurPoint);
//    Point moveKnee = femurHelp.getProjectionPoint(femurKneeCenter);
//    return moveKnee;
//}

Point Knee::getMoveFemurKneeCenter(const FemurImplantInfo& pImplant) const
{
    double resectionThickness;
    Point myInferiorMoveFemurPoint;

    /*
    if (goodSide == MedialSide)
    {
        myInferiorMoveFemurPoint = medialInferiorFemurPoint;
        //resectionThickness = pImplant.femurDistalMedialThickness + LAT_MED_DIFFERENCE - cartilage;
        resectionThickness = pImplant.femurDistalMedialThickness - femurCartilage;
    }
    else
    {
        myInferiorMoveFemurPoint = lateralInferiorFemurPoint;
        //resectionThickness = pImplant.femurDistalLateralThickness - cartilage;
        resectionThickness = pImplant.femurDistalLateralThickness - femurLatMedDiffDistal - femurCartilage;

        //if (resectionThickness < 2)
        //{
         //   resectionThickness = 2;
        //}
    }
    */

    Point femurLineDirectVector = hipCenter - femurKneeCenter;
    femurLineDirectVector.normalice();

    ///////////////////////////////////////////////////////////////////
    Plane helpPlane;
    helpPlane.init(femurLineDirectVector, anteriorCortex);

    if (abs(helpPlane.eval(medialInferiorFemurPoint)) > abs(helpPlane.eval(lateralInferiorFemurPoint)))
    {
        myInferiorMoveFemurPoint = medialInferiorFemurPoint;
        resectionThickness = pImplant.femurDistalMedialThickness - femurCartilage;
    }
    else
    {
        myInferiorMoveFemurPoint = lateralInferiorFemurPoint;
        resectionThickness = pImplant.femurDistalLateralThickness - femurCartilage;
    }
    //////////////////////////////////////////////////////////////////////

    myInferiorMoveFemurPoint = myInferiorMoveFemurPoint + resectionThickness * femurLineDirectVector;

    Plane femurHelp;
    femurHelp.init(femurLineDirectVector, myInferiorMoveFemurPoint);
    Point moveKnee = femurHelp.getProjectionPoint(femurKneeCenter);
    return moveKnee;
}

Point Knee::getPointAtSquareDistance(const Line& line, const Point& pPoint, const Point& nearReferencePoint, float distance, bool closest) const
{
    Point diffPoint = line.getPoint() - pPoint;
    Point equation;
    equation.x = (line.getDirectVector()).dot(line.getDirectVector());
    equation.y = 2.0*((line.getDirectVector()).dot(diffPoint));
    equation.z = (diffPoint.dot(diffPoint)) - (distance*distance);
    double disd = sqrt(equation.y*equation.y - 4.0 * equation.x * equation.z);
    double parameter1 = -(equation.y + disd) / (2.0 * equation.x);
    double parameter2 = -(equation.y - disd) / (2.0 * equation.x);
    Point newPoint1 = (line.getDirectVector() * parameter1) + line.getPoint();
    Point newPoint2 = (line.getDirectVector() * parameter2) + line.getPoint();
    Point diff1 = nearReferencePoint - newPoint1;
    Point diff2 = nearReferencePoint - newPoint2;
    if (closest == ((diff1.x * diff1.x + diff1.y * diff1.y + diff1.z * diff1.z) < (diff2.x * diff2.x + diff2.y * diff2.y + diff2.z * diff2.z)))
        return newPoint1;
    else
        return newPoint2;
}

Point Knee::getPointAtDistance(const Point& directVector, const Point& fixPoint, const Point& nearReferencePoint, float distance, bool closest) const
{
    Point myDirectVector = directVector;
    myDirectVector.normalice();

    Point newPoint1 = fixPoint + (myDirectVector * distance);
    Point newPoint2 = fixPoint - (myDirectVector * distance);
    Point diff1 = nearReferencePoint - newPoint1;
    Point diff2 = nearReferencePoint - newPoint2;
    if (closest == ((diff1.x * diff1.x + diff1.y * diff1.y + diff1.z * diff1.z) < (diff2.x * diff2.x + diff2.y * diff2.y + diff2.z * diff2.z)))
        return newPoint1;
    else
        return newPoint2;
}

void Knee::getDistalInferiorCondyle(Point& lateralPoint, Point& medialPoint)
{
    Point directVectorAxis = hipCenter - femurKneeCenter;
    Point directVectorTEA = lateralEpicondyle - medialEpicondylePerp;
    Point vectorAP = directVectorAxis.cross(directVectorTEA);

    Plane sagitalLateral;
    Plane sagitalMedial;
    Plane transversePlane, transverseCoronal;
    Plane coronal;
    Plane sagital;

    transverseCoronal.init(directVectorAxis, anteriorCortex);
    sagital.init(directVectorTEA, femurKneeCenter);
    coronal.init(vectorAP, medialEpicondyle);
    sagitalLateral.init(directVectorTEA, lateralCondyle);
    sagitalMedial.init(directVectorTEA, medialCondyle);
    transversePlane.init(directVectorAxis, femurKneeCenter + (2.0 * directVectorAxis / sqrt(directVectorAxis.dot(directVectorAxis))));

    if (transversePlane.eval(hipCenter) > 0)
    {
        transversePlane.reverse();
    }

    if (coronal.eval(medialCondyle) > 0)
    {
        coronal.reverse();
    }

    if (sagital.eval(lateralCondyle) < 0)
    {
        sagital.reverse();
    }

    if (transverseCoronal.eval(hipCenter) > 0)
    {
        transverseCoronal.reverse();
    }

    double latDistanceTemp, medDistanceTemp, latDistance, medDistance, coronalLatDistance, coronalMedDistance;

    latDistanceTemp = -1;
    medDistanceTemp = -1;
    coronalLatDistance = -1;
    coronalMedDistance = -1;

    directVectorTEA.normalice();
    Point initLat = lateralCondyle;
    Point initMed = medialCondyle;

    bool fineLat = false, fineMed = false;

    for (double i = -1; i <= 1.01; i += 0.25)
    {
        Point initLatTemp = initLat + i * directVectorTEA;
        Point initMedTemp = initMed + i * directVectorTEA;
        sagitalLateral.movePlane(initLatTemp);
        sagitalMedial.movePlane(initMedTemp);

        vtkSmartPointer<vtkPolyData> contourLat = ImplantTools::getMaxContour(femurPoly, sagitalLateral.getNormalVector(), sagitalLateral.getPoint());
        vtkSmartPointer<vtkPolyData> contourMed = ImplantTools::getMaxContour(femurPoly, sagitalMedial.getNormalVector(), sagitalMedial.getPoint());

        if (!(contourLat->GetPoints()) || !(contourMed->GetPoints()))
        {
            continue;
        }

        if (contourLat->GetPoints()->GetNumberOfPoints() > 0)
        {
            Point pointTemp = ImplantTools::GetFarestPoint(contourLat, transversePlane);
            latDistance = transversePlane.eval(pointTemp);

            if (latDistance > latDistanceTemp)
            {
                latDistanceTemp = latDistance;
                lateralPoint = pointTemp;
                fineLat = true;
            }
        }

        if (contourMed->GetPoints()->GetNumberOfPoints() > 0)
        {
            Point pointTemp = ImplantTools::GetFarestPoint(contourMed, transversePlane);
            medDistance = transversePlane.eval(pointTemp);

            if (medDistance > medDistanceTemp)
            {
                medDistanceTemp = medDistance;
                medialPoint = pointTemp;
                fineMed = true;
            }
        }
    }

    if (fineLat == false || fineMed == false)
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_DISTAL_FEMUR_CONDYLES_POINTS;
    }

    vtkSmartPointer<vtkPoints> allPoints = femurPoly->GetPoints();
    int tSize = allPoints->GetNumberOfPoints();

    for (int i = 0; i < tSize; i++)
    {
        double pnt[3];
        allPoints->GetPoint(i, pnt);
        Point it1(pnt[0], pnt[1], pnt[2]);

        /*if (transversePlane.eval(it1)*transverseSign > 0)
        {
            if (sagitalLateral.isPointNearToPlane(it1, 0.2))
            {
                projP = sagitalLateral.getProjectionPoint(it1);

                latDistance = transversePlane.getDistanceFromPoint(projP);
                if (latDistance > latDistanceTemp)
                {
                    latDistanceTemp = latDistance;
                    lateralPoint = projP;
                }
            }
            else if (sagitalMedial.isPointNearToPlane(it1, 0.2))
            {
                projP = sagitalMedial.getProjectionPoint(it1);

                medDistance = transversePlane.getDistanceFromPoint(projP);
                if (medDistance > medDistanceTemp)
                {
                    medDistanceTemp = medDistance;
                    medialPoint = projP;
                }
            }
        }*/
        if (coronal.eval(it1) > 0 && transverseCoronal.eval(it1) > 0)
        {
            if (sagital.eval(it1) > 2)
            {
                latDistance = coronal.eval(it1);

                if (latDistance > coronalLatDistance)
                {
                    coronalLatDistance = latDistance;
                    coronalDistalLat = (it1);
                }
            }
            else if (sagital.eval(it1) < -2)
            {
                medDistance = coronal.eval(it1);
                if (medDistance > coronalMedDistance)
                {
                    coronalMedDistance = medDistance;
                    coronalDistalMed = (it1);
                }
            }

        }
    }

    Plane lateralPlane, medialPlane, cutPlane;
    medialPlane.init(directVectorTEA, coronalDistalMed);
    lateralPlane.init(directVectorTEA, coronalDistalLat);
    cutPlane = coronal.getPerpendicularPlane(coronalDistalLat, coronalDistalMed);

    if (lateralPlane.eval(medialEpicondylePerp) < 0)
    {
        lateralPlane.reverse();
    }

    if (medialPlane.eval(lateralEpicondyle) < 0)
    {
        medialPlane.reverse();
    }

    vtkSmartPointer<vtkPolyData> contour = ImplantTools::getMaxContour(femurPoly, cutPlane.getNormalVector(), cutPlane.getPoint());

    ////////////////////////////////////////////////////////////////

    Line tempAP(vectorAP, femurKneeCenter);
    Point projPosterior = tempAP.getProjectPoint(medialCondyle);

    vectorAP = femurKneeCenter - projPosterior;
    vectorAP.normalice();

    vtkNew<vtkImplicitPolyDataDistance> polyDistance;
    polyDistance->SetInput(femurPoly);

    Point refPoint = (coronalDistalLat + coronalDistalMed) / 2.;
    Point extPoint = refPoint + 10. * vectorAP;
    Point midGroovePoint;

    ImplantTools::GetInterceptionWithLine(polyDistance, refPoint, extPoint, midGroovePoint);

    vtkIdType refMid = ImplantTools::GetNearestPoints(contour, midGroovePoint);
    vtkIdType refLat = ImplantTools::GetNearestPoints(contour, coronalDistalLat);
    vtkIdType refMed = ImplantTools::GetNearestPoints(contour, coronalDistalMed);

    std::list<std::pair<vtkIdType, vtkIdType>> lines;
    ImplantTools::ExtractSortLines(contour, lines);

    int cont = lines.size();

    while (cont > -1 && lines.front().first != refMid)
    {
        lines.push_back(lines.front());
        lines.pop_front();
        cont--;
    }

    if (cont < 0)
    {
        throw ImplantExceptionCode::CHECK_LANDMARKS_CAN_NOT_DETERMINE_BEGIN_POINT_OF_KNEE_GROOVE;
    }

    std::vector<Point> curve;

    auto it1 = lines.begin();
    auto it2 = lines.end();

    for (; it1 != it2; ++it1)
    {
        if (it1->first != refLat && it1->first != refMed)
        {
            double pnt[3];
            contour->GetPoint((it1)->first, pnt);
            curve.push_back(Point(pnt[0], pnt[1], pnt[2]));
        }
        else
        {
            break;
        }
    }

    auto rit1 = lines.rbegin();
    auto rit2 = lines.rend();

    for (; rit1 != rit2; ++rit1)
    {
        if (rit1->first != refLat && rit1->first != refMed)
        {
            double pnt[3];
            contour->GetPoint((rit1)->first, pnt);
            curve.push_back(Point(pnt[0], pnt[1], pnt[2]));
        }
        else
        {
            break;
        }
    }


    ///////////////////////////////////////////////////////////////

    /*std::list<std::pair<vtkIdType, vtkIdType>> lines;
    ImplantTools::ExtractSortLines(contour, lines);
    auto it1 = lines.begin();
    auto it2 = lines.end();

    std::vector<Point> curve;

    for (; it1 != it2; ++it1)
    {
        double pnt[3];
        contour->GetPoint((it1)->first, pnt);
        if (coronal.eval(pnt) >= 0 && lateralPlane.eval(pnt) > 0 && medialPlane.eval(pnt) > 0)
        {
            curve.push_back(Point(pnt[0], pnt[1], pnt[2]));
        }
    }*/

    Plane basePlane = transversePlane.getPerpendicularPlane(lateralCondyle, medialCondyle);
    basePlane.reverseByPoint(femurKneeCenter);

    Point topPointTemp;

    try
    {
        topPointTemp = ImplantTools::getLocalMinimum(curve, cutPlane, cutPlane.getProjectionVector(basePlane.getNormalVector()));
    }
    catch (...)
    {
        throw ImplantExceptionCode::CHECK_LANDMARKS_CAN_NOT_DETERMINE_BEGIN_POINT_OF_KNEE_GROOVE;
    }

    double myClosest[3];
    double pnt[3] = { topPointTemp.x, topPointTemp.y, topPointTemp.z };

    polyDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);
    topPointOnPatellaPath = Point(myClosest[0], myClosest[1], myClosest[2]);
}

/*
void Knee::getDistalInferiorCondyleOld(Point& lateralPoint, Point& medialPoint) const
{
    Point directVectorAxis = hipCenter - femurKneeCenter;
    Point directVectorTEA = lateralEpicondyle - medialEpicondylePerp;
    Point vectorAP = directVectorAxis.cross(directVectorTEA);

    Plane sagitalLateral;
    Plane sagitalMedial;
    Plane transversePlane;

    sagitalLateral.init(directVectorTEA, lateralCondyle);
    sagitalMedial.init(directVectorTEA, medialCondyle);
    transversePlane.init(directVectorAxis, femurKneeCenter + (2.0 * directVectorAxis / sqrt(directVectorAxis.dot(directVectorAxis))));

    int transverseSign;

    if (transversePlane.eval(hipCenter) > 0)
    {
        transverseSign = -1;
    }
    else
    {
        transverseSign = 1;
    }

    cv::Mat normalPlaneMat(3, 1, CV_64F);
    normalPlaneMat.at <double>(0, 0) = directVectorTEA.x;
    normalPlaneMat.at <double>(1, 0) = directVectorTEA.y;
    normalPlaneMat.at <double>(2, 0) = directVectorTEA.z;

    Point normalXY(0.0, 0.0, 1.0);
    Point rotationAxis = normalXY.cross(directVectorTEA);
    rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
    double rotationAngle = Line::getAngleBetweenVectors(normalXY, directVectorTEA);
    cv::Mat rotation_1 = Line::getRotateMatrix(rotationAxis, -rotationAngle);
    cv::Mat rotation_2 = Line::getRotateMatrix(rotationAxis, rotationAngle);

    cv::Mat rotateVector_1 = rotation_1 * normalPlaneMat;
    cv::Mat rotateVector_2 = rotation_2 * normalPlaneMat;

    cv::Mat rotate;

    double distance_1 = Line::getAngleBetweenVectors(Point(rotateVector_1), normalXY);
    double distance_2 = Line::getAngleBetweenVectors(Point(rotateVector_2), normalXY);

    if (distance_1 < distance_2)
    {
        rotate = rotation_1;
    }
    else
    {
        rotate = rotation_2;
    }

    std::vector<cv::Point2f> coplanarLateral2d, coplanarMedial2d;
    std::vector<Point> lateral3d, medial3d;
    double axisLateralZ = 0.0;
    double axisMedialZ = 0.0;
    auto it1 = mFemur.begin();
    auto it2 = mFemur.end();
    Point projP, rotatePoint;
    cv::Mat rotatePointMat;

    while (it1 != it2)
    {
        if (transversePlane.eval(*it1)*transverseSign > 0)
        {
            if (sagitalLateral.isPointNearToPlane(*it1, 0.2))
            {
                projP = sagitalLateral.getProjectionPoint(*it1);
                lateral3d.push_back(projP);
                rotatePointMat = rotate * projP.ToMatPoint();
                rotatePoint = Point(rotatePointMat);
                coplanarLateral2d.push_back(cv::Point2f(float(rotatePoint.x), float(rotatePoint.y)));
                axisLateralZ += rotatePoint.z;
            }
            else if (sagitalMedial.isPointNearToPlane(*it1, 0.2))
            {
                projP = sagitalMedial.getProjectionPoint(*it1);
                medial3d.push_back(projP);
                rotatePointMat = rotate * projP.ToMatPoint();
                rotatePoint = Point(rotatePointMat);
                coplanarMedial2d.push_back(cv::Point2f(float(rotatePoint.x), float(rotatePoint.y)));
                axisMedialZ += rotatePoint.z;
            }
        }

        ++it1;
    }

    if (coplanarLateral2d.size() * coplanarMedial2d.size() == 0)
    {
        throw ImplantsException("There are no matches between landmarks and the bone surface.");
    }

    axisLateralZ = axisLateralZ / (double(coplanarLateral2d.size()));
    axisMedialZ = axisMedialZ / (double(coplanarMedial2d.size()));

    //cv::RotatedRect boxlateral = cv::minAreaRect(coplanarLateral2d);
    //cv::RotatedRect boxMedial = cv::minAreaRect(coplanarMedial2d);

    cv::RotatedRect boxlateral = cv::fitEllipse(coplanarLateral2d);
    cv::RotatedRect boxMedial = cv::fitEllipse(coplanarMedial2d);

    Line lineAP(vectorAP, femurKneeCenter);
    Point condyleProj = lineAP.getProjectPoint(lateralCondyle);
    Point rightAP = condyleProj - femurKneeCenter;
    rightAP.normalice();
    Point tempLatPoint, tempMedPoint;

    tempLatPoint = getDistalInferiorCondylePoint(boxlateral, rotate, (lateralEpicondyle + 2.0 * rightAP), axisLateralZ);
    tempMedPoint = getDistalInferiorCondylePoint(boxMedial, rotate, (medialEpicondylePerp + 3.0 * rightAP), axisMedialZ);
    directVectorAxis.normalice();

    tempLatPoint = tempLatPoint - 2.0 * directVectorAxis;
    tempMedPoint = tempMedPoint - 2.0 * directVectorAxis;

    it1 = lateral3d.begin();
    it2 = lateral3d.end();

    double maxDistance = Line::getDistanceBetweenPoints(*it1, tempLatPoint, true);
    double tempDistance;
    for (; it1 != it2; ++it1)
    {
        tempDistance = Line::getDistanceBetweenPoints(*it1, tempLatPoint, true);
        if (tempDistance < maxDistance)
        {
            maxDistance = tempDistance;
            lateralPoint = *it1;
        }
    }

    it1 = medial3d.begin();
    it2 = medial3d.end();

    maxDistance = Line::getDistanceBetweenPoints(*it1, tempMedPoint, true);
    for (; it1 != it2; ++it1)
    {
        tempDistance = Line::getDistanceBetweenPoints(*it1, tempMedPoint, true);
        if (tempDistance < maxDistance)
        {
            maxDistance = tempDistance;
            medialPoint = *it1;
        }
    }
}


Point Knee::getDistalInferiorCondylePoint(const cv::RotatedRect& box, const cv::Mat& rotate, const Point& epicondyle, const double axisZ) const
{
    Point directVectorAxis = hipCenter - femurKneeCenter;
    cv::Point2f vertices[4];
    box.points(vertices);

    cv::Mat Va = rotate.inv() * (Point(vertices[0].x, vertices[0].y, axisZ)).ToMatPoint();
    cv::Mat Vb = rotate.inv() * (Point(vertices[1].x, vertices[1].y, axisZ)).ToMatPoint();
    cv::Mat Vc = rotate.inv() * (Point(vertices[2].x, vertices[2].y, axisZ)).ToMatPoint();
    cv::Mat Vd = rotate.inv() * (Point(vertices[3].x, vertices[3].y, axisZ)).ToMatPoint();

    cv::Mat vectorMat1 = Va - Vb;
    cv::Mat vectorMat2 = Vb - Vc;

    Point squareSide1 = Point(vectorMat1);
    Point squareSide2 = Point(vectorMat2);

    double angle1 = (Line::getAngleBetweenVectors(squareSide1, directVectorAxis))*180.0 / PI;
    double angle2 = (Line::getAngleBetweenVectors(squareSide2, directVectorAxis))*180.0 / PI;

    cv::Point3d point1;
    cv::Point3d point2;

    if (abs(90.0 - abs(angle1)) < 50.0 && abs(90.0 - abs(angle1)) > 40.0)
    {
        throw ImplantsException("The lowest condyle cannot be determined automatically.");
    }
    else
    {
        if (abs(90.0 - abs(angle1)) < abs(90.0 - abs(angle2)))
        {
            point1 = Line::getProjectPoint(cv::Point3d(Va), cv::Point3d(Vb), epicondyle);
            point2 = Line::getProjectPoint(cv::Point3d(Vc), cv::Point3d(Vd), epicondyle);
        }
        else
        {
            point1 = Line::getProjectPoint(cv::Point3d(Vb), cv::Point3d(Vc), epicondyle);
            point2 = Line::getProjectPoint(cv::Point3d(Vd), cv::Point3d(Va), epicondyle);
        }

        if (Line::getDistanceBetweenPoints(point1, hipCenter) > Line::getDistanceBetweenPoints(point2, hipCenter))
        {
            return point1;
        }
        else
        {
            return point2;
        }
    }
}

void Knee::get2DSlice(const ImplantImageType::Pointer& imageIn, ImplantImageType::Pointer& imageOut, const int sliceNumber) const
{
    typedef ImplantImageType InputImageType;
    typedef ImplantImageType SliceImageType;
    typedef itk::ExtractImageFilter< InputImageType, SliceImageType > FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetDirectionCollapseToSubmatrix();
    InputImageType::RegionType inputRegion = imageIn->GetLargestPossibleRegion();
    InputImageType::SizeType size = inputRegion.GetSize();
    size[2] = 1;
    InputImageType::IndexType start = inputRegion.GetIndex();
    start[2] = sliceNumber;

    InputImageType::RegionType desiredRegion;
    desiredRegion.SetSize(size);
    desiredRegion.SetIndex(start);
    filter->SetExtractionRegion(desiredRegion);
    filter->SetInput(imageIn);
    try
    {
        filter->Update();
    }
    catch (const itk::ExceptionObject & ex)
    {
        std::cout << ex.GetDescription() << std::endl;
    }
    imageOut = filter->GetOutput();
}*/

std::vector<Point> Knee::getFemurPoints() const
{
    return mFemur;
}
std::vector<Point> Knee::getTibiaPoints() const
{
    return mTibia;
}


Point Knee::getComputeAnkleCenter(const Point& lateralMalleolus, const Point& medialMalleolus)
{
    float k = 44.0 / 100.0; //44% from medial malleolus
    return (medialMalleolus + k * (lateralMalleolus - medialMalleolus));
}

double Knee::getEstimateFemurImplantsDimensions() const
{
    Plane transverse, coronal;
    Point vectorAxis = hipCenter - femurKneeCenter;
    Point vectorML = lateralEpicondyle - medialEpicondylePerp;
    Point vectorAP = vectorAxis.cross(vectorML);
    double distance;
    double implantAverage = 5.0;

    coronal.init(femurDirectVectorAP, femurKneeCenter);
    coronal.movePlane(lateralCondyle);

    if (goodSide == LateralSide)
    {
        transverse.init(vectorAxis, lateralCondyle);
        Point cortexProy = transverse.getProjectionPoint(anteriorCortex);
        Point vectorAPProy = transverse.getProjectionVector(vectorAP);
        Line lateralAP(vectorAPProy, lateralCondyle);
        Point lastCortexProy = lateralAP.getProjectPoint(cortexProy);
        distance = Line::getDistanceBetweenPoints(lastCortexProy, lateralCondyle);// -implantAverage;
    }
    else
    {
        transverse.init(vectorAxis, medialCondyle);
        Point cortexProy = transverse.getProjectionPoint(anteriorCortex);
        Point vectorAPProy = transverse.getProjectionVector(vectorAP);
        Line medialAP(vectorAPProy, medialCondyle);
        Point lastCortexProy = medialAP.getProjectPoint(cortexProy);
        distance = Line::getDistanceBetweenPoints(lastCortexProy, medialCondyle) - abs(coronal.eval(medialCondyle)); //-implantAverage
    }

    return distance;
}

TibiaImplantDimensions Knee::getEstimateTibiaImplantsDimensions() const
{
    //Point freePoint = getTibiaCenterPointFreeAP();
    Plane tibiaCutPlane = getEquisPlaneTibia();
    //tibiaCutPlane.init(tibiaNormalPlaneVector, freePoint);
    //tibiaCutPlane.normalizeNormalVector();

    Point normalXY(0.0, 0.0, 1.0);
    Point rotationAxis = normalXY.cross(tibiaCutPlane.getNormalVector());
    rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
    double rotationAngle = Line::getAngleBetweenVectors(normalXY, tibiaCutPlane.getNormalVector());
    cv::Mat rotation_1 = Line::getRotateMatrix(rotationAxis, -rotationAngle);
    cv::Mat rotation_2 = Line::getRotateMatrix(rotationAxis, rotationAngle);

    cv::Mat rotateVector_1 = rotation_1 * tibiaCutPlane.getNormalVectorMat();
    cv::Mat rotateVector_2 = rotation_2 * tibiaCutPlane.getNormalVectorMat();

    cv::Mat rotate;

    double distance_1 = Line::getAngleBetweenVectors(Point(rotateVector_1), normalXY);
    double distance_2 = Line::getAngleBetweenVectors(Point(rotateVector_2), normalXY);

    if (distance_1 < distance_2)
    {
        rotate = rotation_1;
    }
    else
    {
        rotate = rotation_2;
    }
    std::vector<cv::Point2f> coplanar2d;

    vtkSmartPointer<vtkPolyData> contour = ImplantTools::getContours(tibiaPoly, tibiaCutPlane.getNormalVector(), tibiaCutPlane.getPoint());
    vtkSmartPointer<vtkPoints> contourPoints = contour->GetPoints();
    vtkIdType tSize = contourPoints->GetNumberOfPoints();

    if (tSize < 25)
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_DIMENSIONS_OF_TIBIA_IMPLANT;
    }

    double x, y;

    for (vtkIdType i = 0; i < tSize; i++)
    {
        double pnt[3];
        contourPoints->GetPoint(i, pnt);

        Point pointProj(pnt[0], pnt[1], pnt[2]);
        cv::Mat point(3, 1, CV_64F);
        point.at <double>(0, 0) = pointProj.x;
        point.at <double>(1, 0) = pointProj.y;
        point.at <double>(2, 0) = pointProj.z;
        cv::Mat rotatePointMat = rotate * point;
        x = rotatePointMat.at <double>(0, 0);
        y = rotatePointMat.at <double>(1, 0);
        coplanar2d.push_back(cv::Point2f(float(x), float(y)));
    }

    /*auto it1 = mTibia.begin();
    auto it2 = mTibia.end();
    double x, y;

    for (; it1 != it2; ++it1)
    {
        if (tibiaCutPlane.isPointNearToPlane(*it1, 1.0) == true)
        {
            Point pointProj = tibiaCutPlane.getProjectionPoint(*it1);
            cv::Mat point(3, 1, CV_64F);
            point.at <double>(0, 0) = pointProj.x;
            point.at <double>(1, 0) = pointProj.y;
            point.at <double>(2, 0) = pointProj.z;
            cv::Mat rotatePointMat = rotate * point;
            x = rotatePointMat.at <double>(0, 0);
            y = rotatePointMat.at <double>(1, 0);
            coplanar2d.push_back(cv::Point2f(float(x), float(y)));
        }
    }*/

    //cv::RotatedRect box = cv::fitEllipse(coplanar2d);

    cv::RotatedRect box = cv::minAreaRect(coplanar2d);

    double with = box.size.width;
    double height = box.size.height;
    TibiaImplantDimensions result;

    if (with > height)
    {
        result.long_axis = with;
        result.short_axis = height;
    }
    else
    {
        result.long_axis = height;
        result.short_axis = with;
    }

    return result;
}

Side Knee::getGoodSide(const Point& hipCenter, const Point& kneeCenter, const Point& latEpi, const Point& medEpi, const Point& ankle) const
{

    Point femurAxis = hipCenter - kneeCenter;
    Point vectorTEA = latEpi - medEpi;
    Plane sagital, transverse;
    transverse.init(femurAxis, kneeCenter);
    sagital.init(vectorTEA, kneeCenter);
    Point ankleProj = transverse.getProjectionPoint(ankle);

    if (sagital.eval(latEpi) < 0)
    {
        sagital.reverse();
    }

    double eval = sagital.eval(ankleProj);

    if (eval > 0)
    {
        return MedialSide;
    }
    else
    {
        return LateralSide;
    }

    /*Point femurAxis = hipCenter - kneeCenter;
    Point legAxis = hipCenter - ankle;
    Plane transverse, coronal;
    transverse.init(femurAxis, kneeCenter);
    coronal = transverse.getPerpendicularPlane(latEpi, medEpi);

    Point hipProj = coronal.getProjectionPoint(hipCenter);
    Point latProj = coronal.getProjectionPoint(latEpi);
    Point medProj = coronal.getProjectionPoint(medEpi);
    Point kneeProj = coronal.getProjectionPoint(kneeCenter);
    Point legAxisProj = coronal.getProjectionVector(legAxis);

    Line lineRef = Line(legAxisProj, hipProj);
    Point epiVector = latProj - medProj;
    epiVector.normalice();

    Point latRef = kneeProj + epiVector;
    Point medRef = kneeProj - epiVector;

    if (lineRef.getDistanceFromPoint(latRef) > lineRef.getDistanceFromPoint(medRef))
    {
        return LateralSide;
    }
    else
    {
        return MedialSide;
    }*/
}

vtkSmartPointer<vtkPolyData> Knee::GetFemurPoly() const
{
    return femurPoly;
}

vtkSmartPointer<vtkPolyData> Knee::GetTibiaPoly() const
{
    return tibiaPoly;
}

vtkSmartPointer<vtkPolyData> Knee::GetPatellaPoly() const
{
    return mPatella.getPatellaPoly();
}

double Knee::getFemurCartilage() const
{
    return femurCartilage;
}

double Knee::getTibiaCartilage() const
{
    return tibiaCartilage;
}

bool Knee::getIsVarus() const
{
    return isVarus;
}

//Point Knee::getLateralCoronalDistalPoint() const
//{
//    return coronalDistalLat;
//}
//
//Point Knee::getMedialCoronalDistalPoint() const
//{
//    return coronalDistalMed;
//}
//
//Point Knee::getTopPointOnPatellaPath() const
//{
//    return topPointOnPatellaPath;
//}

void Knee::getAutomaticPlateaus()
{
    Point axis = tibiaKneeCenter - ankleCenter;
    Plane axial;
    axial.init(axis, tibiaTubercle);

    Point vectorAP = tibiaTubercle - axial.getProjectionPoint(pclCenter);
    axis.normalice();
    vectorAP.normalice();
    Point vectorTEA;

    TemplateTibia templateObj;
    Point sourceAP, sourceTEA, sourceAxis, sourcePoint;

    std::vector<cv::Point3d> myTemplatePoints;

    Point latPlateauIn, latPlateauOut, medPlateauIn, medPlateauOut;

    if (isRight == true)
    {
        vectorTEA = vectorAP.cross(axis);

        sourceAP = templateObj.vectorRightAP;
        sourceTEA = templateObj.vectorRightTEA;
        sourcePoint = templateObj.kneeCenterRight;

        myTemplatePoints = templateObj.mTemplateRight;

        latPlateauIn = templateObj.rightLatPlateauIn;
        latPlateauOut = templateObj.rightLatPlateauOut;
        medPlateauIn = templateObj.rightMedPlateauIn;
        medPlateauOut = templateObj.rightMedPlateauOut;
    }
    else
    {
        vectorTEA = axis.cross(vectorAP);

        sourceAP = templateObj.vectorLeftAP;
        sourceTEA = templateObj.vectorLeftTEA;
        sourcePoint = templateObj.kneeCenterLeft;

        myTemplatePoints = templateObj.mTemplateLeft;

        latPlateauIn = templateObj.leftLatPlateauIn;
        latPlateauOut = templateObj.leftLatPlateauOut;
        medPlateauIn = templateObj.leftMedPlateauIn;
        medPlateauOut = templateObj.leftMedPlateauOut;
    }
    vectorTEA.normalice();

    sourceAxis = sourceAP.cross(sourceTEA);
    sourceAxis.normalice();

    Point targetAxis = vectorAP.cross(vectorTEA);
    targetAxis.normalice();

    std::vector<cv::Point3d> vectorTarget = { vectorAP, vectorTEA, targetAxis };
    std::vector<cv::Point3d> vectorSource = { sourceAP, sourceTEA, sourceAxis };

    cv::Mat data(7, 1, CV_64F);
    double tScale = 1.2;

    cv::Mat rotation = LeastSquaresScaleICP::GetRotationAnglesXYZ(vectorSource, vectorTarget, data);
    cv::Mat translation = tibiaKneeCenter.ToMatPoint() - (rotation * sourcePoint.ToMatPoint());

    data.at<double>(3, 0) = data.at<double>(0, 0);
    data.at<double>(4, 0) = data.at<double>(1, 0);
    data.at<double>(5, 0) = data.at<double>(2, 0);

    data.at<double>(0, 0) = translation.at<double>(0, 0);
    data.at<double>(1, 0) = translation.at<double>(1, 0);
    data.at<double>(2, 0) = translation.at<double>(2, 0);

    data.at<double>(6, 0) = tScale;

    LeastSquaresScaleICP registerObj(myTemplatePoints);

    double error = registerObj.LeastSquaresScale(tibiaPoly, data);

    cv::Mat myTranslation(3, 1, CV_64F);
    myTranslation.at<double>(0, 0) = data.at<double>(0, 0);
    myTranslation.at<double>(1, 0) = data.at<double>(1, 0);
    myTranslation.at<double>(2, 0) = data.at<double>(2, 0);

    double angleX = data.at<double>(3, 0);
    double angleY = data.at<double>(4, 0);
    double angleZ = data.at<double>(5, 0);

    double scale = data.at<double>(6, 0);

    cv::Mat myRotation = registerObj.GetRotationMatrix(angleX, angleY, angleZ);

    cv::Mat transformPointMatIn = scale * (myRotation * latPlateauIn.ToMatPoint()) + myTranslation;
    cv::Mat transformPointMatOut = scale * (myRotation * latPlateauOut.ToMatPoint()) + myTranslation;
    
    latPlateauIn = Point(transformPointMatIn);
    latPlateauOut = Point(transformPointMatOut);

    transformPointMatIn = scale * (myRotation * medPlateauIn.ToMatPoint()) + myTranslation;
    transformPointMatOut = scale * (myRotation * medPlateauOut.ToMatPoint()) + myTranslation;

    medPlateauIn = Point(transformPointMatIn);
    medPlateauOut = Point(transformPointMatOut);

    vtkNew<vtkImplicitPolyDataDistance> polyDistance;
    polyDistance->SetInput(tibiaPoly);

    double dist1 = ImplantTools::GetInterceptionWithLine(polyDistance, latPlateauIn, latPlateauOut, this->lateralPlateau);
    double dist2 = ImplantTools::GetInterceptionWithLine(polyDistance, medPlateauIn, medPlateauOut, this->medialPlateau);

    double pntLat[3] = { this->lateralPlateau.x, this->lateralPlateau.y, this->lateralPlateau.z };
    double pntMed[3] = { this->medialPlateau.x, this->medialPlateau.y, this->medialPlateau.z };

    double myClosestLat[3];
    double myClosestMed[3];

    polyDistance->EvaluateFunctionAndGetClosestPoint(pntLat, myClosestLat);
    polyDistance->EvaluateFunctionAndGetClosestPoint(pntMed, myClosestMed);

    this->lateralPlateau.x = myClosestLat[0];
    this->lateralPlateau.y = myClosestLat[1];
    this->lateralPlateau.z = myClosestLat[2];

    this->medialPlateau.x = myClosestMed[0];
    this->medialPlateau.y = myClosestMed[1];
    this->medialPlateau.z = myClosestMed[2];
}

void Knee::makeKneeGroovePath()
{
    Point vectorAxis = hipCenter - femurKneeCenter;
    Point vectorTEA = lateralEpicondyle - medialEpicondylePerp;
    Point vectorAP = vectorAxis.cross(vectorTEA);

    Point distalCenter = (coronalDistalLat + coronalDistalMed) / 2.0;

    std::vector<Point> points;
    Plane coronal, sagital, axial;

    coronal.init(vectorAP, femurKneeCenter);
    coronal.reverseByPoint(medialCondyle, false);
    axial.init(vectorAxis, coronalDistalLat);
    axial.reverseByPoint(hipCenter, false);

    sagital = coronal.getPerpendicularPlane(femurKneeCenter, distalCenter);

    vtkSmartPointer<vtkPolyData> contour = ImplantTools::getMaxContour(femurPoly, sagital.getNormalVector(), sagital.getPoint());

    if (contour->GetNumberOfPoints() == 0)
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_PATELLA_GROOVE;
    }

    std::list<std::pair<vtkIdType, vtkIdType>> lines;
    ImplantTools::ExtractSortLines(contour, lines);
    auto it1 = lines.begin();
    auto it2 = lines.end();
    for (; it1 != it2; ++it1)
    {
        double pnt[3];
        contour->GetPoint((it1)->first, pnt);
        if (coronal.eval(pnt) >= 0 && axial.eval(pnt) > 0)
        {
            points.push_back(Point(pnt[0], pnt[1], pnt[2]));
        }
    }

    ////////////////////////////////////////////////////////////////////////

    Plane basePlane = axial.getPerpendicularPlane(lateralCondyle, medialCondyle);
    basePlane.reverseByPoint(femurKneeCenter);

    std::vector<Point> resultTemp;

    Plane latPlane, medPlane;
    latPlane = axial.getPerpendicularPlane(lateralCondyle, coronalDistalLat);
    medPlane = axial.getPerpendicularPlane(medialCondyle, coronalDistalMed);
    latPlane.reverseByPoint(medialEpicondylePerp);
    medPlane.reverseByPoint(lateralEpicondyle);

    int step = points.size() / 30;
    if (step == 0)
    {
        step = 1;
    }

    Plane axialTemp;
    Point axialTempPoint = femurKneeCenter + 0.3 * (hipCenter - femurKneeCenter);

    axialTemp.init(vectorAxis, axialTempPoint);
    Point pointFit1 = axialTemp.getProjectionPoint(medialCondyle);
    Point pointFit2 = axialTemp.getProjectionPoint(lateralCondyle);

    int contTemp = 0;
    for (int i = step; i < points.size(); i += step)
    {
        Plane tempPlane;
        std::list<Point> parabola;
        tempPlane.init(pointFit1, pointFit2, points[i]);

        vtkSmartPointer<vtkPolyData> contour = ImplantTools::getMaxContour(femurPoly, tempPlane.getNormalVector(), tempPlane.getPoint());

        if (contour->GetNumberOfPoints() == 0)
        {
            continue;
        }

        vtkIdType pointRef = ImplantTools::GetNearestPoints(contour, points[i]);

        std::list<std::pair<vtkIdType, vtkIdType>> lines;
        ImplantTools::ExtractSortLines(contour, lines);

        int cont = lines.size();

        while (cont > -1 && lines.front().first != pointRef)
        {
            lines.push_back(lines.front());
            lines.pop_front();
            cont--;
        }

        if (cont < 0)
        {
            continue;
        }

        auto it1 = lines.begin();
        auto it2 = lines.end();

        for (; it1 != it2; ++it1)
        {
            double pnt[3];
            contour->GetPoint((it1)->first, pnt);

            if (latPlane.eval(pnt) > 0 && medPlane.eval(pnt) > 0)
            {
                parabola.push_back(Point(pnt[0], pnt[1], pnt[2]));
            }
            else
            {
                break;
            }
        }

        auto rit1 = lines.rbegin();
        auto rit2 = lines.rend();

        for (; rit1 != rit2; ++rit1)
        {
            double pnt[3];
            contour->GetPoint((rit1)->first, pnt);

            if (latPlane.eval(pnt) > 0 && medPlane.eval(pnt) > 0)
            {
                parabola.push_front(Point(pnt[0], pnt[1], pnt[2]));
            }
            else
            {
                break;
            }
        }

        Point myVertex;

        try
        {
            myVertex = ImplantTools::getLocalMinimum(parabola, tempPlane, tempPlane.getProjectionVector(basePlane.getNormalVector()));
        }
        catch (...)
        {
            throw ImplantExceptionCode::CHECK_LANDMARKS_CAN_NOT_DETERMINE_WHITESIDE_LINE;
        }

        resultTemp.push_back(myVertex);
    }

    bool resultOk;
    Plane fixPlane = Plane::getBestPlaneOutliers(resultTemp, mKneeGrooveOutLiers, resultOk, 0.4);

    if (resultOk == false)
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_PATELLA_GROOVE;
    }

    fixPlane.reverseByPoint(medialEpicondyle);
    Point condyleVector = medialCondyle - lateralCondyle;
    Plane helpPlane;
    helpPlane.init(vectorAxis, femurKneeCenter);
    condyleVector = helpPlane.getProjectionVector(condyleVector);
    Point vectorEpi = medialEpicondylePerp - lateralEpicondyle;
    Point vectorGroove = fixPlane.getNormalVector();
    vectorGroove = helpPlane.getProjectionVector(vectorGroove);

    double angleEpi = ImplantTools::getAngleBetweenVectorsDegree(condyleVector, vectorEpi);
    double angleGroove = ImplantTools::getAngleBetweenVectorsDegree(condyleVector, vectorGroove);

    // Fixing medialEpicondylePerp attribute
    if (abs(angleGroove - 3.0) < abs(angleEpi - 3.0))
    {
        Line Epi(vectorGroove, lateralEpicondyle);
        this->medialEpicondylePerp = Epi.getProjectPoint(medialEpicondyle);
    }

    vtkSmartPointer<vtkPolyData> lastContour = ImplantTools::getMaxContour(femurPoly, fixPlane.getNormalVector(), fixPlane.getPoint());
    mKneeGroovePath.clear();
    std::list<std::pair<vtkIdType, vtkIdType>> lastLines;
    ImplantTools::ExtractSortLines(lastContour, lastLines);
    it1 = lastLines.begin();
    it2 = lastLines.end();
    for (; it1 != it2; ++it1)
    {
        double pnt[3];
        lastContour->GetPoint((it1)->first, pnt);
        if (coronal.eval(pnt) >= 0 && axial.eval(pnt) > 0)
        {
            mKneeGroovePath.push_back(Point(pnt[0], pnt[1], pnt[2]));
        }
    }
}

Plane Knee::getPatellaFrontPlane() const
{
    return mPatellaPlane;
}

Point Knee::getPatellaCenter() const
{
    return mPatella.getPatellaCenter();
}

Point Knee::getPatellaDistalPosteriorPoint() const
{
    return mPatellaDistalPosteriorPoint;
}

std::vector<Point> Knee::getKneeGroovePath() const
{
    return mKneeGroovePath;
}

//std::vector<Point> Knee::getKneeGrooveOutLiers() const
//{
//    return mKneeGrooveOutLiers;
//}

double Knee::getPatellaThickness() const
{
    return mPatella.getPatellaThickness();
}

double Knee::getPatellaDiameter() const
{
    return mPatella.getPatellaDiameter();
}

Patella Knee::getPatella() const
{
    return mPatella;
}

void Knee::setFemurCartilage(double pCartilage)
{
    femurCartilage = pCartilage;
}

void Knee::setTibiaCartilage(double pCartilage)
{
    tibiaCartilage = pCartilage;
}

double Knee::getInitialAnglePCA(const Point& hipCenter, const Point& femurKneeCenter, const Point& lateralEpicondyle,
    const Point& medialEpicondyle, const Point& lateralCondyle, const Point& medialCondyle)
{
    Point forceLine = hipCenter - femurKneeCenter;
    Point teaLine = medialEpicondyle - lateralEpicondyle;
    Point condyleLine = medialCondyle - lateralCondyle;

    Plane axial;
    axial.init(forceLine, femurKneeCenter);

    Point vectorTEA = axial.getProjectionVector(teaLine);
    Point vectorCondyle = axial.getProjectionVector(condyleLine);

    return ImplantTools::getAngleBetweenVectorsDegree(vectorTEA, vectorCondyle);
}



