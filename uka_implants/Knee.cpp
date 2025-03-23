
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
#include "LeastSquaresScaleICP.hpp"
#include "TemplatePoints.hpp"

using namespace UKA::IMPLANTS;

Knee::Knee()
{
    isInit = false;
}

void Knee::init(const Point& hipCenter, const Point& femurKneeCenter, const Point& lateralEpicondyle,
    const Point& medialEpicondyle, const Point& tibiaKneeCenter, const Point& tibiaTubercle, const Point& pclCenter, 
    const Point& ankleCenter, const vtkSmartPointer<vtkPolyData> femurPoly, 
    const vtkSmartPointer<vtkPolyData> tibiaPoly, KneeSideEnum pSide, SurgerySideEnum pSurgery, bool findRefPoints, double cartilage, uint8_t imageValueMax)
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
    this->tibiaKneeCenter = tibiaKneeCenter;
    this->tibiaTubercle = tibiaTubercle;
	this->mSurgerySize = pSurgery;

    Point forceLineVector = hipCenter - femurKneeCenter;
	forceLineVector.normalice();
    Point TEA = lateralEpicondyle - medialEpicondylePerp;
	TEA.normalice();

	this->cortexRef = femurKneeCenter + (ImplantTools::getDistanceBetweenPoints(lateralEpicondyle, medialEpicondyle)/2) * forceLineVector;

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

	if (findRefPoints == true)
	{
		findAutomaticPlateaus();
		findFemurCondylesUsingLeastSquare();
	}

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
    isInit = true;
}


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

void Knee::findFemurCondylesUsingDistalPoints()
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

    transverse.init(forceLine, cortexRef);
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

SurgerySideEnum Knee::getSurgerySide() const
{
	return mSurgerySize;
}

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

Point Knee::getTibiaVectorToSurgicalSideTEA() const
{
	if (mSurgerySize == SurgerySideEnum::KLateral)
	{
		return getTibiaVectorLateralTEA();
	}
	else
	{
		return -getTibiaVectorLateralTEA();
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

Point Knee::getCortexRef() const
{
    return cortexRef;
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

    if (mSurgerySize == SurgerySideEnum::KMedial)
    {
        myMoveCondyle = medialCondyle;
    }
    else
    {
        myMoveCondyle = lateralCondyle;
    }
    //////////////////////////////////////////////////////////////////////
	//std::cout << "--------------------------->>>>>> " << resectionThickness << std::endl;
	//std::cout << "--------------------------->>>>>> Thickness " << pImplant.femurPosteriorThickness << " Cartilago: "<< femurCartilage << std::endl;
	resectionThickness = pImplant.femurPosteriorThickness - femurCartilage;
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

Point Knee::getMovePlateauCenter(const TibiaImplantInfo& pImplant) const
{
	double resectionThickness;
	Point myMovePlateau;
	
	Point tibiaLineDirectVector = tibiaNormalPlaneVector;
	tibiaLineDirectVector.normalice();

	///////////////////////////////////////////////////////////////////

	Point lateralCenter, medialCenter;
	getTibiaPlateauCenters(lateralCenter, medialCenter);

	if (mSurgerySize == SurgerySideEnum::KMedial)
	{
		myMovePlateau = medialCenter;
	}
	else
	{
		myMovePlateau = lateralCenter;
	}

	//////////////////////////////////////////////////////////////////////
	resectionThickness = pImplant.tibiaThickness - tibiaCartilage;
	myMovePlateau = myMovePlateau - resectionThickness * tibiaLineDirectVector;
	return myMovePlateau;
}

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

    if (mSurgerySize == SurgerySideEnum::KMedial)
    {
        myMovePlateau = medialPlateau;
    }
    else
    {
        myMovePlateau = lateralPlateau;
    }
    //////////////////////////////////////////////////////////////////////
	resectionThickness = pImplant.tibiaThickness - tibiaCartilage;
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

    /*Point forceLine = hipCenter - femurKneeCenter;
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
    }*/

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

    if (mSurgerySize == SurgerySideEnum::KMedial)
    {
        myInferiorMoveFemurPoint = medialInferiorFemurPoint;
    }
    else
    {
        myInferiorMoveFemurPoint = lateralInferiorFemurPoint;
    }
    //////////////////////////////////////////////////////////////////////
	resectionThickness = pImplant.femurDistalThickness - femurCartilage;
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

    if (mSurgerySize == SurgerySideEnum::KMedial)
    {
        myMovePlateau = medialPlateau;
        resectionThickness = pImplant.tibiaThickness - tibiaCartilage;
    }
    else
    {
        myMovePlateau = lateralPlateau;
        resectionThickness = pImplant.tibiaThickness - tibiaCartilage;
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

void Knee::getTibiaPlateauCenters(Point& lateral, Point& medial) const
{
	Plane equis = getEquisPlaneTibia();

	vtkSmartPointer<vtkPolyData> contourMax = ImplantTools::getContours(GetTibiaPoly(), equis.getNormalVector(), equis.getPoint());
	vtkSmartPointer<vtkPoints> vtkMyPoints = contourMax->GetPoints();
	int tVtkPointSize = vtkMyPoints->GetNumberOfPoints();

	std::vector<Point> contourLat, cotourMed;
	Plane midPlane = equis.getPerpendicularPlane(equis.getProjectionPoint(pclCenter), equis.getProjectionPoint(tibiaTubercle));
	Point tempPoint = midPlane.getPoint();
	midPlane.movePlane(medialPlateau);
	midPlane.reverseByPoint(lateralPlateau);
	midPlane.movePlane(tempPoint);

	Point refLateralVector = lateralPlateau - equis.getProjectionPoint(lateralPlateau);
	Point refMedialVector = medialPlateau - equis.getProjectionPoint(medialPlateau);

	for (int i = 0; i < tVtkPointSize; i++)
	{
		double pnt[3];
		vtkMyPoints->GetPoint(i, pnt);
		Point currentPoint(pnt[0], pnt[1], pnt[2]);
		if (midPlane.eval(currentPoint) >= 1)
		{
			contourLat.push_back(currentPoint + refLateralVector);
		}
		else if (midPlane.eval(currentPoint) <= -2)
		{
			cotourMed.push_back(currentPoint + refMedialVector);
		}
	}

	if (contourLat.size() > 0 && cotourMed.size() > 0)
	{
		lateral = ImplantTools::getPolygonCenter(contourLat, equis.getNormalVector());
		medial = ImplantTools::getPolygonCenter(cotourMed, equis.getNormalVector());
	}
	else
	{
		lateral = lateralPlateau;
		medial = medialPlateau;
	}
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

    if (mSurgerySize == SurgerySideEnum::KMedial)
    {
        myInferiorMoveFemurPoint = medialInferiorFemurPoint;
        resectionThickness = pImplant.femurDistalThickness - femurCartilage;
    }
    else
    {
        myInferiorMoveFemurPoint = lateralInferiorFemurPoint;
        resectionThickness = pImplant.femurDistalThickness - femurCartilage;
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

    transverseCoronal.init(directVectorAxis, cortexRef);
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
}

Point Knee::getComputeAnkleCenter(const Point& lateralMalleolus, const Point& medialMalleolus)
{
    float k = 44.0 / 100.0; //44% from medial malleolus
    return (medialMalleolus + k * (lateralMalleolus - medialMalleolus));
}

double Knee::getEstimateFemurImplantsDimensions() const
{
 //   Plane transverse, coronal;
 //   Point vectorAxis = hipCenter - femurKneeCenter;
 //   Point vectorML = lateralEpicondyle - medialEpicondylePerp;
 //   Point vectorAP = vectorAxis.cross(vectorML);
 //   double distance;
 //   //double implantAverage = 5.0;

 //   coronal.init(femurDirectVectorAP, femurKneeCenter);
 //   coronal.movePlane(lateralCondyle);
	//double epicondyleRadiusDistance = ImplantTools::getDistanceBetweenPoints(lateralEpicondyle, lateralCondyle);

    if (goodSide == LateralSide)
    {
        //transverse.init(vectorAxis, lateralCondyle);
        //Point cortexProy = transverse.getProjectionPoint(anteriorCortex);
        //Point vectorAPProy = transverse.getProjectionVector(vectorAP);
        //Line lateralAP(vectorAPProy, lateralCondyle);
        //Point lastCortexProy = lateralAP.getProjectPoint(cortexProy);
        //distance = Line::getDistanceBetweenPoints(lastCortexProy, lateralCondyle);// -implantAverage;
		return 2 * ImplantTools::getDistanceBetweenPoints(lateralEpicondyle, lateralCondyle);
    }
    else
    {
        //transverse.init(vectorAxis, medialCondyle);
        //Point cortexProy = transverse.getProjectionPoint(anteriorCortex);
        //Point vectorAPProy = transverse.getProjectionVector(vectorAP);
        //Line medialAP(vectorAPProy, medialCondyle);
        //Point lastCortexProy = medialAP.getProjectPoint(cortexProy);
        //distance = Line::getDistanceBetweenPoints(lastCortexProy, medialCondyle) - abs(coronal.eval(medialCondyle)); //-implantAverage
		
		return 2 * ImplantTools::getDistanceBetweenPoints(medialEpicondyle, medialCondyle);
	}

    //return distance;
}

/*
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
*/

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
}

vtkSmartPointer<vtkPolyData> Knee::GetFemurPoly() const
{
    return femurPoly;
}

vtkSmartPointer<vtkPolyData> Knee::GetTibiaPoly() const
{
    return tibiaPoly;
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

void Knee::findAutomaticPlateaus()
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

void Knee::findFemurCondylesUsingLeastSquare()
{
	vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
	implicitPolyDataDistance->SetInput(GetFemurPoly());

	TemplateFemur templateObj;

	Point targetAP, targetTEA, targetAxis, targetCenter;
	Point sourceAP, sourceTEA, sourceAxis, sourceCenter;
	double tScale;

	targetAP = getFemurDirectVectorAP();
	targetTEA = getFemurVectorLateralTEA();
	targetAxis = targetAP.cross(targetTEA);
	targetAxis.normalice();
	targetCenter = (getLateralEpicondyle() + getMedialEpicondyle()) / 2.0;

	Point diff = getLateralEpicondyle() - getMedialEpicondyle();
	double targetSize = sqrt(diff.dot(diff));

	std::vector<cv::Point3d> myTemplatePoints;

	Point latInferiorCondyleIn, medInferiorCondyleIn, latPosteriorCondyleIn, medPosteriorCondyleIn;
	Point latInferiorCondyleOut, medInferiorCondyleOut, latPosteriorCondyleOut, medPosteriorCondyleOut;

	if (getIsRight() == false)
	{
		sourceAP = templateObj.vectorLeftAP;
		sourceTEA = templateObj.vectorLeftTEA;
		sourceCenter = templateObj.centerLeft;
		tScale = (targetSize / templateObj.sizeLeft);
		myTemplatePoints = templateObj.mTemplateLeft;

		latInferiorCondyleIn = templateObj.leftLatInferiorCondyleIn;
		medInferiorCondyleIn = templateObj.leftMedInferiorCondyleIn;
		latPosteriorCondyleIn = templateObj.leftLatPosteriorCondyleIn;
		medPosteriorCondyleIn = templateObj.leftMedPosteriorCondyleIn;

		latInferiorCondyleOut = templateObj.leftLatInferiorCondyleOut;
		medInferiorCondyleOut = templateObj.leftMedInferiorCondyleOut;
		latPosteriorCondyleOut = templateObj.leftLatPosteriorCondyleOut;
		medPosteriorCondyleOut = templateObj.leftMedPosteriorCondyleOut;
	}
	else
	{
		sourceAP = templateObj.vectorRightAP;
		sourceTEA = templateObj.vectorRightTEA;
		sourceCenter = templateObj.centerRight;
		tScale = (targetSize / templateObj.sizeLeft);
		myTemplatePoints = templateObj.mTemplateRight;

		latInferiorCondyleIn = templateObj.rightLatInferiorCondyleIn;
		medInferiorCondyleIn = templateObj.rightMedInferiorCondyleIn;
		latPosteriorCondyleIn = templateObj.rightLatPosteriorCondyleIn;
		medPosteriorCondyleIn = templateObj.rightMedPosteriorCondyleIn;

		latInferiorCondyleOut = templateObj.rightLatInferiorCondyleOut;
		medInferiorCondyleOut = templateObj.rightMedInferiorCondyleOut;
		latPosteriorCondyleOut = templateObj.rightLatPosteriorCondyleOut;
		medPosteriorCondyleOut = templateObj.rightMedPosteriorCondyleOut;
	}

	if (tScale < 1.0)
	{
		tScale = 1.0;
	}

	sourceAxis = sourceAP.cross(sourceTEA);
	sourceAxis.normalice();

	std::vector<cv::Point3d> vectorTarget = { targetAP, targetTEA, targetAxis };
	std::vector<cv::Point3d> vectorSource = { sourceAP, sourceTEA, sourceAxis };

	cv::Mat data(7, 1, CV_64F);

	cv::Mat rotation = LeastSquaresScaleICP::GetRotationAnglesXYZ(vectorSource, vectorTarget, data);
	cv::Mat translation = targetCenter.ToMatPoint() - (rotation * sourceCenter.ToMatPoint());

	data.at<double>(3, 0) = data.at<double>(0, 0);
	data.at<double>(4, 0) = data.at<double>(1, 0);
	data.at<double>(5, 0) = data.at<double>(2, 0);

	data.at<double>(0, 0) = translation.at<double>(0, 0);
	data.at<double>(1, 0) = translation.at<double>(1, 0);
	data.at<double>(2, 0) = translation.at<double>(2, 0);

	data.at<double>(6, 0) = tScale;

	LeastSquaresScaleICP registerObj(myTemplatePoints);

	registerObj.LeastSquaresScale(GetFemurPoly(), data);

	////////////////////////////////////////////////////////////////////

	cv::Mat myTranslation(3, 1, CV_64F);
	myTranslation.at<double>(0, 0) = data.at<double>(0, 0);
	myTranslation.at<double>(1, 0) = data.at<double>(1, 0);
	myTranslation.at<double>(2, 0) = data.at<double>(2, 0);

	double angleX = data.at<double>(3, 0);
	double angleY = data.at<double>(4, 0);
	double angleZ = data.at<double>(5, 0);

	double scale = data.at<double>(6, 0);

	cv::Mat myRotation = registerObj.GetRotationMatrix(angleX, angleY, angleZ);

	std::vector<cv::Mat> pointsIn, pointsOut;

	pointsIn.push_back(scale * (myRotation * latInferiorCondyleIn.ToMatPoint()) + myTranslation);
	pointsIn.push_back(scale * (myRotation * medInferiorCondyleIn.ToMatPoint()) + myTranslation);
	pointsIn.push_back(scale * (myRotation * latPosteriorCondyleIn.ToMatPoint()) + myTranslation);
	pointsIn.push_back(scale * (myRotation * medPosteriorCondyleIn.ToMatPoint()) + myTranslation);

	pointsOut.push_back(scale * (myRotation * latInferiorCondyleOut.ToMatPoint()) + myTranslation);
	pointsOut.push_back(scale * (myRotation * medInferiorCondyleOut.ToMatPoint()) + myTranslation);
	pointsOut.push_back(scale * (myRotation * latPosteriorCondyleOut.ToMatPoint()) + myTranslation);
	pointsOut.push_back(scale * (myRotation * medPosteriorCondyleOut.ToMatPoint()) + myTranslation);

	for (int i = 0; i < 4; i++)
	{
		Point temp;
		ImplantTools::GetInterceptionWithLine(implicitPolyDataDistance, pointsIn[i], pointsOut[i], temp);
		double pnt[3] = { temp.x, temp.y, temp.z };
		double closest[3];
		implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, closest);

		if (i == 0)
		{
			lateralInferiorFemurPoint.x = closest[0];
			lateralInferiorFemurPoint.y = closest[1];
			lateralInferiorFemurPoint.z = closest[2];
		}
		else if (i == 1)
		{
			medialInferiorFemurPoint.x = closest[0];
			medialInferiorFemurPoint.y = closest[1];
			medialInferiorFemurPoint.z = closest[2];
		}
		else if (i == 2)
		{
			lateralCondyle.x = closest[0];
			lateralCondyle.y = closest[1];
			lateralCondyle.z = closest[2];
		}
		else
		{
			medialCondyle.x = closest[0];
			medialCondyle.y = closest[1];
			medialCondyle.z = closest[2];
		}
	}

}



