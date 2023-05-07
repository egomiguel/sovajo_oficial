
#include "vtkImplicitPolyDataDistance.h"
#include "vtkClipClosedSurface.h"
#include "vtkPlaneCollection.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkMatrix4x4.h"
#include "vtkClipPolyData.h"
#include "vtkFeatureEdges.h"
#include "vtkStripper.h"
#include "vtkAppendPolyData.h"
#include "Line.hpp"
#include "Balance.hpp"
#include <itkMatrix.h>
#include "ImplantsException.hpp"
#include "vtkCellLocator.h"
#include "ImplantTools.hpp"

using namespace TKA::IMPLANTS;

double balanceSquareDistance(const Point& a, const Point& b)
{
    Point diff = a - b;
    return diff.dot(diff);
}

Balance::Balance()
{
    isInit = false;
}

void Balance::init(const Knee& pKnee, const FemurImplant& pFemurImplant, const TibiaImplant& pTibiaImplant, const itk::Rigid3DTransform<>::Pointer pImplantToBoneFemurTransform, const itk::Rigid3DTransform<>::Pointer pImplantToBoneTibiaTransform)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_BALANCE;
    }
    this->knee_ = pKnee;
    this->femurTransformCtToMarker = cv::Mat::eye(4, 4, CV_64FC1);
    this->femurTransformMarkerToCamera = cv::Mat::eye(4, 4, CV_64FC1);
    this->tibiaTransformCtToMarker = cv::Mat::eye(4, 4, CV_64FC1);
    this->tibiaTransformMarkerToCamera = cv::Mat::eye(4, 4, CV_64FC1);
    this->kneeCapTransformCtToMarker = cv::Mat::eye(4, 4, CV_64FC1);
    this->kneeCapTransformMarkerToCamera = cv::Mat::eye(4, 4, CV_64FC1);
    this->femurImplant = pFemurImplant;
    this->tibiaImplant = pTibiaImplant;

    mFemurTransformImplantToBone = Rigid3DTransformToCV(pImplantToBoneFemurTransform);
    mTibiaTransformImplantToBone = Rigid3DTransformToCV(pImplantToBoneTibiaTransform);

    femurImplantPoly = pFemurImplant.GetTransformImplantModel(pImplantToBoneFemurTransform);

    double resizeVector = 1000000.0;
    Point anteriorPoint, posteriorPoint;

    anteriorPoint = knee_.getFemurKneeCenter() + resizeVector * knee_.getFemurDirectVectorAP();
    posteriorPoint = knee_.getFemurKneeCenter() - resizeVector * knee_.getFemurDirectVectorAP();

    PlaneA = TransformImplantPlaneToBone(pFemurImplant.getPlaneA(), pImplantToBoneFemurTransform);
    PlaneA.reverseByPoint(anteriorPoint);

    PlaneB = TransformImplantPlaneToBone(pFemurImplant.getPlaneB(), pImplantToBoneFemurTransform);
    PlaneB.reverseByPoint(anteriorPoint);

    PlaneC = TransformImplantPlaneToBone(pFemurImplant.getPlaneC(), pImplantToBoneFemurTransform);
    PlaneC.reverseByPoint(knee_.getHipCenter());

    PlaneD = TransformImplantPlaneToBone(pFemurImplant.getPlaneD(), pImplantToBoneFemurTransform);
    PlaneD.reverseByPoint(posteriorPoint);

	PlaneMid = TransformImplantPlaneToBone(pFemurImplant.getMidPlane(), pImplantToBoneFemurTransform);

    PlaneTibia = TransformImplantPlaneToBone(pTibiaImplant.getTibiaPlane(), pImplantToBoneTibiaTransform);
    PlaneTibia.reverseByPoint(knee_.getAnkleCenter(), false);

    //setImplantPlateaus(tibiaImplant.getFirstPlateau(), tibiaImplant.getSecondPlateau(), PlaneTibia, pImplantToBoneTibiaTransform);

    cv::Point3d planeNormal, planePoint;

    vtkNew<vtkPlane> vtkPlaneA, vtkPlaneB, vtkPlaneC, vtkPlaneD, vtkPlaneTibia;

    planeNormal = PlaneA.getNormalVector();
    planePoint = PlaneA.getPoint();
    vtkPlaneA->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
    vtkPlaneA->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);

    planeNormal = PlaneB.getNormalVector();
    planePoint = PlaneB.getPoint();
    vtkPlaneB->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
    vtkPlaneB->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);

    planeNormal = PlaneC.getNormalVector();
    planePoint = PlaneC.getPoint();
    vtkPlaneC->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
    vtkPlaneC->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);

    planeNormal = PlaneD.getNormalVector();
    planePoint = PlaneD.getPoint();
    vtkPlaneD->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
    vtkPlaneD->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);

    vtkNew<vtkPlaneCollection> femurPlanes;
    femurPlanes->AddItem(vtkPlaneA);
    femurPlanes->AddItem(vtkPlaneB);
    femurPlanes->AddItem(vtkPlaneC);
    femurPlanes->AddItem(vtkPlaneD);

    vtkNew<vtkClipClosedSurface> femurClipper;
    femurClipper->SetInputData(knee_.GetFemurPoly());
    femurClipper->SetClippingPlanes(femurPlanes);
    femurClipper->Update();

    femurPolyCut = femurClipper->GetOutput();

    //femurCoordenate = new CoordenateSystemFemur(pKnee.getHipCenter(), pKnee.getLateralEpicondyle(), pKnee.getMedialEpicondylePerp(), pKnee.getFemurKneeCenter());

    isInit = true;
}

vtkSmartPointer<vtkPolyData> Balance::getFemurPolyCut() const
{
    return femurPolyCut;
}

Balance::~Balance()
{
   /* delete femurCoordenate;
    femurCoordenate = NULL;*/
}

vtkSmartPointer<vtkPolyData> Balance::TransformPolyFemurToCamera(const vtkSmartPointer<vtkPolyData> poly) const
{
    cv::Mat transform = (femurTransformMarkerToCamera * femurTransformCtToMarker);

    vtkNew<vtkMatrix4x4> m;
    m->SetElement(0, 0, transform.at<double>(0, 0));
    m->SetElement(0, 1, transform.at<double>(0, 1));
    m->SetElement(0, 2, transform.at<double>(0, 2));
    m->SetElement(0, 3, transform.at<double>(0, 3));
    m->SetElement(1, 0, transform.at<double>(1, 0));
    m->SetElement(1, 1, transform.at<double>(1, 1));
    m->SetElement(1, 2, transform.at<double>(1, 2));
    m->SetElement(1, 3, transform.at<double>(1, 3));
    m->SetElement(2, 0, transform.at<double>(2, 0));
    m->SetElement(2, 1, transform.at<double>(2, 1));
    m->SetElement(2, 2, transform.at<double>(2, 2));
    m->SetElement(2, 3, transform.at<double>(2, 3));
    m->SetElement(3, 0, transform.at<double>(3, 0));
    m->SetElement(3, 1, transform.at<double>(3, 1));
    m->SetElement(3, 2, transform.at<double>(3, 2));
    m->SetElement(3, 3, transform.at<double>(3, 3));

    vtkNew<vtkTransform> vtkTransform;
    vtkTransform->SetMatrix(m);

    vtkNew<vtkTransformPolyDataFilter> transformFilter;
    transformFilter->SetInputData(poly);
    transformFilter->SetTransform(vtkTransform);
    transformFilter->Update();
    return transformFilter->GetOutput();
}

vtkSmartPointer<vtkPolyData> Balance::TransformPolyTibiaToCamera(const vtkSmartPointer<vtkPolyData> poly) const
{
    cv::Mat transform = (tibiaTransformMarkerToCamera * tibiaTransformCtToMarker);

    vtkNew<vtkMatrix4x4> m;
    m->SetElement(0, 0, transform.at<double>(0, 0));
    m->SetElement(0, 1, transform.at<double>(0, 1));
    m->SetElement(0, 2, transform.at<double>(0, 2));
    m->SetElement(0, 3, transform.at<double>(0, 3));
    m->SetElement(1, 0, transform.at<double>(1, 0));
    m->SetElement(1, 1, transform.at<double>(1, 1));
    m->SetElement(1, 2, transform.at<double>(1, 2));
    m->SetElement(1, 3, transform.at<double>(1, 3));
    m->SetElement(2, 0, transform.at<double>(2, 0));
    m->SetElement(2, 1, transform.at<double>(2, 1));
    m->SetElement(2, 2, transform.at<double>(2, 2));
    m->SetElement(2, 3, transform.at<double>(2, 3));
    m->SetElement(3, 0, transform.at<double>(3, 0));
    m->SetElement(3, 1, transform.at<double>(3, 1));
    m->SetElement(3, 2, transform.at<double>(3, 2));
    m->SetElement(3, 3, transform.at<double>(3, 3));

    vtkNew<vtkTransform> vtkTransform;
    vtkTransform->SetMatrix(m);

    vtkNew<vtkTransformPolyDataFilter> transformFilter;
    transformFilter->SetInputData(poly);
    transformFilter->SetTransform(vtkTransform);
    transformFilter->Update();
    return transformFilter->GetOutput();
}

PointTypeITK Balance::getKneeCapPoint() const
{
    Point kneeCap = knee_.getPatellaRotationPoint();

    cv::Mat mat(4, 1, CV_64FC1);
    mat.at <double>(0, 0) = kneeCap.x;
    mat.at <double>(1, 0) = kneeCap.y;
    mat.at <double>(2, 0) = kneeCap.z;
    mat.at <double>(3, 0) = 1.0;

    cv::Mat resultMat = (kneeCapTransformMarkerToCamera * kneeCapTransformCtToMarker) * mat;
    Point newKneeCapPoint(resultMat.at<double>(0, 0), resultMat.at<double>(1, 0), resultMat.at<double>(2, 0));

    PointTypeITK itkPoint;
    itkPoint[0] = newKneeCapPoint.x;
    itkPoint[1] = newKneeCapPoint.y;
    itkPoint[2] = newKneeCapPoint.z;

    return itkPoint;
}

Point Balance::transformFemurPointToCamera(const Point& phisicalPoint) const
{
    cv::Mat mat(4, 1, CV_64FC1);
    mat.at <double>(0, 0) = phisicalPoint.x;
    mat.at <double>(1, 0) = phisicalPoint.y;
    mat.at <double>(2, 0) = phisicalPoint.z;
    mat.at <double>(3, 0) = 1.0;

    cv::Mat result = (femurTransformMarkerToCamera * femurTransformCtToMarker) * mat;
    return Point(result.at<double>(0, 0), result.at<double>(1, 0), result.at<double>(2, 0));
}

Point Balance::transformFemurVectorToCamera(const Point& pLandmarkVector) const
{
    cv::Mat transformMat = femurTransformMarkerToCamera * femurTransformCtToMarker;

    double* matrix = new double[9];

    matrix[0] = transformMat.at<double>(0, 0);
    matrix[1] = transformMat.at<double>(0, 1);
    matrix[2] = transformMat.at<double>(0, 2);

    matrix[3] = transformMat.at<double>(1, 0);
    matrix[4] = transformMat.at<double>(1, 1);
    matrix[5] = transformMat.at<double>(1, 2);

    matrix[6] = transformMat.at<double>(2, 0);
    matrix[7] = transformMat.at<double>(2, 1);
    matrix[8] = transformMat.at<double>(2, 2);

    cv::Mat rotationMat(3, 3, CV_64FC1, matrix);
    cv::Mat newVector = rotationMat * pLandmarkVector.ToMatPoint();
    return Point(newVector);
}

Point Balance::transformTibiaVectorToCamera(const Point& pLandmarkVector) const
{
    cv::Mat transformMat = tibiaTransformMarkerToCamera * tibiaTransformCtToMarker;

    double* matrix = new double[9];

    matrix[0] = transformMat.at<double>(0, 0);
    matrix[1] = transformMat.at<double>(0, 1);
    matrix[2] = transformMat.at<double>(0, 2);

    matrix[3] = transformMat.at<double>(1, 0);
    matrix[4] = transformMat.at<double>(1, 1);
    matrix[5] = transformMat.at<double>(1, 2);

    matrix[6] = transformMat.at<double>(2, 0);
    matrix[7] = transformMat.at<double>(2, 1);
    matrix[8] = transformMat.at<double>(2, 2);

    cv::Mat rotationMat(3, 3, CV_64FC1, matrix);
    cv::Mat newVector = rotationMat * pLandmarkVector.ToMatPoint();
    return Point(newVector);
}

//Plane Balance::transformFemurPlane(const Plane& femurCT) const
//{
//    Point pPoint = femurCT.getPoint();
//    cv::Mat normalMat = femurCT.getNormalVectorMat();
//
//    cv::Mat pointMat(4, 1, CV_64FC1);
//    pointMat.at <double>(0, 0) = pPoint.x;
//    pointMat.at <double>(1, 0) = pPoint.y;
//    pointMat.at <double>(2, 0) = pPoint.z;
//    pointMat.at <double>(3, 0) = 1.0;
//
//    cv::Mat transformMat = femurTransformMarkerToCamera * femurTransformCtToMarker;
//
//    double* matrix = new double[9];
//
//    matrix[0] = transformMat.at<double>(0, 0);
//    matrix[1] = transformMat.at<double>(0, 1);
//    matrix[2] = transformMat.at<double>(0, 2);
//
//    matrix[3] = transformMat.at<double>(1, 0);
//    matrix[4] = transformMat.at<double>(1, 1);
//    matrix[5] = transformMat.at<double>(1, 2);
//
//    matrix[6] = transformMat.at<double>(2, 0);
//    matrix[7] = transformMat.at<double>(2, 1);
//    matrix[8] = transformMat.at<double>(2, 2);
//
//    cv::Mat rotationMat(3, 3, CV_64FC1, matrix);
//
//    cv::Mat newNormal = rotationMat * normalMat;
//    cv::Mat newPoint = transformMat * pointMat;
//
//    Point finalPoint(newPoint.at<double>(0, 0), newPoint.at<double>(1, 0), newPoint.at<double>(2, 0));
//    Point finalNormal = Point(newNormal);
//
//    Plane finalPlane;
//    finalPlane.init(finalNormal, finalPoint);
//    finalPlane.normalizeNormalVector();
//    return finalPlane;
//}
//
Plane Balance::transformTibiaPlane(const Plane& tibiaCT) const
{
    Point pPoint = tibiaCT.getPoint();
    cv::Mat normalMat = tibiaCT.getNormalVectorMat();

    cv::Mat pointMat(4, 1, CV_64FC1);
    pointMat.at <double>(0, 0) = pPoint.x;
    pointMat.at <double>(1, 0) = pPoint.y;
    pointMat.at <double>(2, 0) = pPoint.z;
    pointMat.at <double>(3, 0) = 1.0;

    cv::Mat transformMat = tibiaTransformMarkerToCamera * tibiaTransformCtToMarker;

    double* matrix = new double[9];

    matrix[0] = transformMat.at<double>(0, 0);
    matrix[1] = transformMat.at<double>(0, 1);
    matrix[2] = transformMat.at<double>(0, 2);

    matrix[3] = transformMat.at<double>(1, 0);
    matrix[4] = transformMat.at<double>(1, 1);
    matrix[5] = transformMat.at<double>(1, 2);

    matrix[6] = transformMat.at<double>(2, 0);
    matrix[7] = transformMat.at<double>(2, 1);
    matrix[8] = transformMat.at<double>(2, 2);

    cv::Mat rotationMat(3, 3, CV_64FC1, matrix);

    cv::Mat newNormal = rotationMat * normalMat;
    cv::Mat newPoint = transformMat * pointMat;

    Point finalPoint(newPoint.at<double>(0, 0), newPoint.at<double>(1, 0), newPoint.at<double>(2, 0));
    Point finalNormal = Point(newNormal);

    Plane finalPlane;
    finalPlane.init(finalNormal, finalPoint);
    finalPlane.normalizeNormalVector();
    return finalPlane;
}

Point Balance::transformTibiaPointToCamera(const Point& phisicalPoint) const
{
    cv::Mat mat(4, 1, CV_64FC1);
    mat.at <double>(0, 0) = phisicalPoint.x;
    mat.at <double>(1, 0) = phisicalPoint.y;
    mat.at <double>(2, 0) = phisicalPoint.z;
    mat.at <double>(3, 0) = 1.0;

    cv::Mat result = (tibiaTransformMarkerToCamera * tibiaTransformCtToMarker) * mat;
    return Point(result.at<double>(0, 0), result.at<double>(1, 0), result.at<double>(2, 0));
}

Plane Balance::TransformImplantPlaneToBone(const Plane& plane, const itk::Rigid3DTransform<>::Pointer transform) const
{
    cv::Mat myRotation = Rigid3DTransformToCVRotation(transform);
    cv::Mat myTranslation = Rigid3DTransformToCVTranslation(transform);

    cv::Mat transformNormalVector = myRotation * plane.getNormalVectorMat();
    cv::Mat transformPoint = (myRotation * plane.getPointMat()) + myTranslation;
    Plane transformPlane;
    transformPlane.init(Point(transformNormalVector), Point(transformPoint));
    return transformPlane;
}

void Balance::setImplantPlateaus(const Point& implantPlateau1, const Point& implantPlateau2, const Plane& tibia, const itk::Rigid3DTransform<>::Pointer transform)
{
    cv::Mat myRotation = Rigid3DTransformToCVRotation(transform);
    cv::Mat myTranslation = Rigid3DTransformToCVTranslation(transform);

    cv::Mat pointResult1 = myRotation * (implantPlateau1.ToMatPoint()) + myTranslation;

    cv::Mat pointResult2 = myRotation * (implantPlateau2.ToMatPoint()) + myTranslation;

    Point implantPlateauT1 = Point(pointResult1);

    Point implantPlateauT2 = Point(pointResult2);

    Point medialProj = tibia.getProjectionPoint(knee_.getMedialPlateau());

}

cv::Mat Balance::Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer transform) const
{
    itk::Matrix< double, 3, 3 > rotation = transform->GetMatrix();
    double* matrix = new double[9];

    matrix[0] = rotation[0][0];
    matrix[1] = rotation[0][1];
    matrix[2] = rotation[0][2];

    matrix[3] = rotation[1][0];
    matrix[4] = rotation[1][1];
    matrix[5] = rotation[1][2];

    matrix[6] = rotation[2][0];
    matrix[7] = rotation[2][1];
    matrix[8] = rotation[2][2];

    cv::Mat matrixCV(3, 3, CV_64FC1, matrix);
    return matrixCV;
}

cv::Mat Balance::Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer transform) const
{
    itk::Vector< double, 3 > translate = transform->GetOffset();
    cv::Mat result(3, 1, CV_64FC1);

    result.at <double>(0, 0) = translate[0];
    result.at <double>(1, 0) = translate[1];
    result.at <double>(2, 0) = translate[2];
    return result;
}

cv::Mat Balance::Rigid3DTransformToCV(const itk::Rigid3DTransform<>::Pointer transform) const
{
    itk::Matrix< double, 3, 3 > rotation = transform->GetMatrix();
    itk::Vector< double, 3 > translate = transform->GetOffset();
    double* matrix = new double[16];

    matrix[0] = rotation[0][0];
    matrix[1] = rotation[0][1];
    matrix[2] = rotation[0][2];
    matrix[3] = translate[0];

    matrix[4] = rotation[1][0];
    matrix[5] = rotation[1][1];
    matrix[6] = rotation[1][2];
    matrix[7] = translate[1];

    matrix[8] = rotation[2][0];
    matrix[9] = rotation[2][1];
    matrix[10] = rotation[2][2];
    matrix[11] = translate[2];

    matrix[12] = 0;
    matrix[13] = 0;
    matrix[14] = 0;
    matrix[15] = 1;

    cv::Mat matrixCV(4, 4, CV_64FC1, matrix);
    return matrixCV;
}

void Balance::setTransformFemurCtToMarker(const itk::Rigid3DTransform<>::Pointer transform)
{
    femurTransformCtToMarker = Rigid3DTransformToCV(transform);
}

void Balance::setTransformTibiaCtToMarker(const itk::Rigid3DTransform<>::Pointer transform)
{
    tibiaTransformCtToMarker = Rigid3DTransformToCV(transform);
}

void Balance::setTransformFemurMarkerToCamera(const itk::Rigid3DTransform<>::Pointer transform)
{
    femurTransformMarkerToCamera = Rigid3DTransformToCV(transform);
}

void Balance::setTransformTibiaMarkerToCamera(const itk::Rigid3DTransform<>::Pointer transform)
{
    tibiaTransformMarkerToCamera = Rigid3DTransformToCV(transform);
}

void Balance::setTransformKneeCapCtToMarker(const itk::Rigid3DTransform<>::Pointer transform)
{
    kneeCapTransformCtToMarker = Rigid3DTransformToCV(transform);
}

void Balance::setTransformKneeCapMarkerToCamera(const itk::Rigid3DTransform<>::Pointer transform)
{
    kneeCapTransformMarkerToCamera = Rigid3DTransformToCV(transform);
}

ImplantImageType::PointType Balance::cvPointToITK(const Point& phisicalPoint) const
{
    ImplantImageType::PointType itkPoint;
    itkPoint[0] = phisicalPoint.x;
    itkPoint[1] = phisicalPoint.y;
    itkPoint[2] = phisicalPoint.z;
    return itkPoint;
}

Point Balance::itkPointToCV(const ImplantImageType::PointType& phisicalPoint) const
{
    Point cvPoint;
    cvPoint.x = phisicalPoint[0];
    cvPoint.y = phisicalPoint[1];
    cvPoint.z = phisicalPoint[2];
    return cvPoint;
}

BalanceInfo Balance::distanceByAngleBeforeResectionBone() const
{
    double distanceLateral;
    double distanceMedial;

    Point closest;
    Point latPlateau = transformTibiaPointToCamera(knee_.getLateralPlateau());
    Point medPlateau = transformTibiaPointToCamera(knee_.getMedialPlateau());

    vtkSmartPointer<vtkPolyData> poly = TransformPolyFemurToCamera(knee_.GetFemurPoly());

    double cartilage = 0;// knee_.getFemurCartilage() + knee_.getTibiaCartilage();

    distanceLateral = closestPoint(poly, latPlateau, closest) - cartilage;
    distanceMedial = closestPoint(poly, medPlateau, closest) - cartilage;

    return BalanceInfo(distanceLateral, distanceMedial);
}

BalanceInfo Balance::distanceByAngleWithImplant(const PointTypeITK& implantLateralPlateau, const PointTypeITK& implantMedialPlateau, bool useImplantThickness) const
{
    double distanceLateral;
    double distanceMedial;

    Point closest;

    Point latPlateau(implantLateralPlateau[0], implantLateralPlateau[1], implantLateralPlateau[2]);
    Point medPlateau(implantMedialPlateau[0], implantMedialPlateau[1], implantMedialPlateau[2]);

    //Point latPlateau = transformTibiaPoint(implantPlateauLat);
    //Point medPlateau = transformTibiaPoint(implantPlateauMed);

    vtkSmartPointer<vtkPolyData> polyImplantFemur = TransformPolyFemurToCamera(femurImplantPoly);

    double thicknessLat = 0;
    double thicknessMed = 0;

    if (useImplantThickness == true)
    {
        thicknessLat = femurImplant.getImplantInfo().femurDistalLateralThickness + tibiaImplant.getImplantInfo().tibiaLateralThickness;
        thicknessMed = femurImplant.getImplantInfo().femurDistalMedialThickness + tibiaImplant.getImplantInfo().tibiaMedialThickness;
    }

    distanceLateral = closestPoint(polyImplantFemur, latPlateau, closest) + thicknessLat;
    distanceMedial = closestPoint(polyImplantFemur, medPlateau, closest) + thicknessMed;

    return BalanceInfo(distanceLateral, distanceMedial);
}

Plane Balance::ComputeNewPlaneTibia(const itk::Rigid3DTransform<>::Pointer pImplantToBoneTibiaTransform)
{
    Plane PlaneTibiaNew = TransformImplantPlaneToBone(tibiaImplant.getTibiaPlane(), pImplantToBoneTibiaTransform);
    return PlaneTibiaNew;
}

vtkSmartPointer<vtkPolyData> Balance::ComputeNewPolyFemur(const itk::Rigid3DTransform<>::Pointer pImplantToBoneFemurTransform)
{
    Point anteriorPoint, posteriorPoint;
    anteriorPoint = knee_.getFemurKneeCenter() + 1000 * knee_.getFemurDirectVectorAP();
    posteriorPoint = knee_.getFemurKneeCenter() - 1000 * knee_.getFemurDirectVectorAP();

    Plane newPlaneA = TransformImplantPlaneToBone(femurImplant.getPlaneA(), pImplantToBoneFemurTransform);
    newPlaneA.reverseByPoint(anteriorPoint);

    Plane newPlaneB = TransformImplantPlaneToBone(femurImplant.getPlaneB(), pImplantToBoneFemurTransform);
    newPlaneB.reverseByPoint(anteriorPoint);

    Plane newPlaneC = TransformImplantPlaneToBone(femurImplant.getPlaneC(), pImplantToBoneFemurTransform);
    newPlaneC.reverseByPoint(knee_.getHipCenter());

    Plane newPlaneD = TransformImplantPlaneToBone(femurImplant.getPlaneD(), pImplantToBoneFemurTransform);
    newPlaneD.reverseByPoint(posteriorPoint);

    cv::Point3d planeNormal, planePoint;

    vtkNew<vtkPlane> vtkPlaneA, vtkPlaneB, vtkPlaneC, vtkPlaneD, vtkPlaneTibia;

    planeNormal = newPlaneA.getNormalVector();
    planePoint = newPlaneA.getPoint();
    vtkPlaneA->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
    vtkPlaneA->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);

    planeNormal = newPlaneB.getNormalVector();
    planePoint = newPlaneB.getPoint();
    vtkPlaneB->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
    vtkPlaneB->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);

    planeNormal = newPlaneC.getNormalVector();
    planePoint = newPlaneC.getPoint();
    vtkPlaneC->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
    vtkPlaneC->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);

    planeNormal = newPlaneD.getNormalVector();
    planePoint = newPlaneD.getPoint();
    vtkPlaneD->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
    vtkPlaneD->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);

    vtkNew<vtkPlaneCollection> femurPlanes;
    femurPlanes->AddItem(vtkPlaneA);
    femurPlanes->AddItem(vtkPlaneB);
    femurPlanes->AddItem(vtkPlaneC);
    femurPlanes->AddItem(vtkPlaneD);

    vtkNew<vtkClipClosedSurface> femurClipper;
    femurClipper->SetInputData(knee_.GetFemurPoly());
    femurClipper->SetClippingPlanes(femurPlanes);
    femurClipper->Update();

    femurPolyCut = femurClipper->GetOutput();

    return femurPolyCut;
}

BalanceInfo Balance::distanceByAngleFemurImplantToTibiaPlane() const
{
	Point midPlaneVector = transformFemurVectorToCamera(PlaneMid.getNormalVector());
	Point midPlanePoint = transformFemurPointToCamera(PlaneMid.getPoint());
	Plane implantMidPlane;
	implantMidPlane.init(midPlaneVector, midPlanePoint);

	Point latEpicondyle = transformFemurPointToCamera(knee_.getLateralEpicondyle());
	Point medEpicondyle = transformFemurPointToCamera(knee_.getMedialEpicondyle());

	vtkSmartPointer<vtkPolyData> poly = TransformPolyFemurToCamera(femurImplantPoly);

	implantMidPlane.reverseByPoint(latEpicondyle);

	auto latData = ImplantTools::GetDistancePlaneToSurface(poly, PlaneTibia, implantMidPlane);
	implantMidPlane.reverse();
	auto medData = ImplantTools::GetDistancePlaneToSurface(poly, PlaneTibia, implantMidPlane);

	return BalanceInfo(latData.first, medData.first);
}

BalanceInfo Balance::distanceByAngleAfterResectionBone(const Plane& pTibia, const vtkSmartPointer<vtkPolyData> pFemurPoly, double toolSize) const
{
    double distanceLateral;
    double distanceMedial;

    Point closest;

    Point latPlateau = pTibia.getProjectionPoint(knee_.getLateralPlateau());
    Point medPlateau = pTibia.getProjectionPoint(knee_.getMedialPlateau());

    latPlateau = transformTibiaPointToCamera(latPlateau);
    medPlateau = transformTibiaPointToCamera(medPlateau);

    vtkSmartPointer<vtkPolyData> poly = TransformPolyFemurToCamera(pFemurPoly);

    distanceLateral = closestPoint(poly, latPlateau, closest) - toolSize;
    distanceMedial = closestPoint(poly, medPlateau, closest) - toolSize;

    return BalanceInfo(distanceLateral, distanceMedial);
}

BalanceInfo Balance::distanceByAngleAfterResectionBone(double toolSize) const
{
    double distanceLateral;
    double distanceMedial;

    Point closest;

    Point latPlateau = PlaneTibia.getProjectionPoint(knee_.getLateralPlateau());
    Point medPlateau = PlaneTibia.getProjectionPoint(knee_.getMedialPlateau());

    latPlateau = transformTibiaPointToCamera(latPlateau);
    medPlateau = transformTibiaPointToCamera(medPlateau);

    vtkSmartPointer<vtkPolyData> poly = TransformPolyFemurToCamera(femurPolyCut);

    distanceLateral = closestPoint(poly, latPlateau, closest) - toolSize;
    distanceMedial = closestPoint(poly, medPlateau, closest) - toolSize;

    return BalanceInfo(distanceLateral, distanceMedial);
}

AnglesInfo Balance::anglesByMotion() const
{
    LegAngle leg;
    Point ankleCenter = transformTibiaPointToCamera(knee_.getAnkleCenter());
    Point tibiaKneeCenter = transformTibiaPointToCamera(knee_.getTibiaKneeCenter());
    Point femurKneeCenter = transformFemurPointToCamera(knee_.getFemurKneeCenter());
    Point hipCenter = transformFemurPointToCamera(knee_.getHipCenter());
    Point medialCondyle = transformFemurPointToCamera(knee_.getMedialCondyle());

    Point lateral = transformFemurPointToCamera(knee_.getLateralEpicondyle());
    Point medial = transformFemurPointToCamera(knee_.getMedialEpicondylePerp());

    Point femurAP = transformFemurVectorToCamera(knee_.getFemurDirectVectorAP());
    Point femurTEA = transformFemurVectorToCamera(knee_.getFemurVectorLateralTEA());
    Point tibiaTEA = transformTibiaVectorToCamera(knee_.getTibiaVectorLateralTEA());
	Point tibiaAP = transformTibiaVectorToCamera(knee_.getTibiaDirectVectorAP());////////////////////////

    double flexion_angle = leg.getFlexionAngle(hipCenter, femurKneeCenter, lateral, medial, tibiaKneeCenter, ankleCenter, knee_.getIsRight());
    double varus_angle = leg.getVarusAngle(femurKneeCenter, medialCondyle, femurAP, femurTEA, tibiaTEA, tibiaAP, knee_.getIsRight());
    double rotation_angle = leg.getRotationAngle(hipCenter, femurKneeCenter, femurTEA, tibiaTEA, tibiaAP, ankleCenter, knee_.getIsRight());

    return AnglesInfo(flexion_angle, varus_angle, rotation_angle);
}

vtkSmartPointer<vtkPolyData> Balance::GetFemurPolyCut() const
{
    return femurPolyCut;
}

double Balance::closestPoint(const vtkSmartPointer<vtkPolyData> poly, const Point& p, Point& closest) const
{
    double point[3] = { p.x, p.y, p.z };
    return closestPoint(poly, point, closest);
}

double Balance::closestPoint(const vtkSmartPointer<vtkPolyData> poly, double * point, Point& closest) const
{
    /*vtkNew<vtkCellLocator> cellLocator;
    cellLocator->SetDataSet(poly);
    cellLocator->BuildLocator();

    double closestPoint[3];
    double squareDistance;
    vtkIdType cellId;
    int subId;

    cellLocator->FindClosestPoint(point, closestPoint, cellId, subId, squareDistance);

    closest = VtkToCv(closestPoint);

    return sqrt(squareDistance);*/

    vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
    implicitPolyDataDistance->SetInput(poly);
    double myClosest[3];
    double signedDistance = implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(point, myClosest);
    closest = Point(myClosest[0], myClosest[1], myClosest[2]);
    return signedDistance;
}

double Balance::GetTriangleArea(const Point& a, const Point& b, const Point& c) const
{
    Point v1 = b - a;
    Point v2 = c - a;
    Point normal = v1.cross(v2);
    double area = sqrt(normal.dot(normal));
    return area / 2.0;
}

double Balance::GetClosestPointOnTriangle(const Point& a, const Point& b, const Point& c, const Point& myPoint, Point& closestPoint) const
{
    Plane plane;
    plane.init(a, b, c);
    Point projectPoint = plane.getProjectionPoint(myPoint);

    double bigArea = GetTriangleArea(a, b, c);
    double totalArea = GetTriangleArea(a, b, projectPoint) + GetTriangleArea(b, c, projectPoint) + GetTriangleArea(c, a, projectPoint);

    if (abs(totalArea - bigArea) <= EPSILON)
    {
        closestPoint = projectPoint;
        return Line::getDistanceBetweenPoints(projectPoint, myPoint, true);
    }
    double nearDistance, distance;
    Point nearPoint, nextP1, nextP2;

    distance = Line::getDistanceBetweenPoints(a, projectPoint, true);
    nearDistance = distance;
    nearPoint = a;
    nextP1 = b;
    nextP2 = c;

    distance = Line::getDistanceBetweenPoints(b, projectPoint, true);

    if (distance < nearDistance)
    {
        nearDistance = distance;
        nearPoint = b;
        nextP1 = a;
        nextP2 = c;
    }

    distance = Line::getDistanceBetweenPoints(c, projectPoint, true);

    if (distance < nearDistance)
    {
        nearPoint = c;
        nextP1 = a;
        nextP2 = b;
    }
    Point closestPoint1, closestPoint2;
    double distance1 = GetDistancePointToSegment(nearPoint, nextP1, projectPoint, closestPoint1);
    double distance2 = GetDistancePointToSegment(nearPoint, nextP2, projectPoint, closestPoint2);

    if (distance1 < distance2)
    {
        closestPoint = closestPoint1;
        return Line::getDistanceBetweenPoints(closestPoint1, myPoint, true);
    }
    else
    {
        closestPoint = closestPoint2;
        return Line::getDistanceBetweenPoints(closestPoint2, myPoint, true);
    }
}

double Balance::GetDistancePointToSegment(const Point& S1, const Point& S2, const Point& P, Point& closestPoint) const
{
    Line nearLine = Line::makeLineWithPoints(S1, S2);
    Point projOnLine = nearLine.getProjectPoint(P);

    Point vectorRef = S2 - S1;
    Point Ref2 = S2 - projOnLine;
    Point Ref1 = projOnLine - S1;

    if (Ref2.dot(vectorRef) < 0)
    {
        closestPoint = S2;
    }
    else if (Ref1.dot(vectorRef) < 0)
    {
        closestPoint = S1;
    }
    else
    {
        closestPoint = projOnLine;
    }
    return Line::getDistanceBetweenPoints(closestPoint, P, true);
}

cv::Point3d Balance::VtkToCv(double * P) const
{
    cv::Point3d result = cv::Point3d(P[0], P[1], P[2]);
    return result;
}

itk::Rigid3DTransform<>::Pointer Balance::getNewImplantToBoneFemurTransformAxial(const std::vector<PointTypeITK>& axialPointsCT)
{
    bool result;
    Plane axial = Plane::getBestPlane(axialPointsCT, result);
    axial.reverseByPoint(knee_.getHipCenter());

    Point newNormal = axial.getNormalVector();
    Point oldNormal = PlaneC.getNormalVector();
    Point oldPoint = PlaneC.getPoint();

    cv::Mat lastRotation = ImplantTools::GetGeneralRotateTransformVectors(oldNormal, newNormal);
    cv::Mat lastTranslation = axial.getPoint().ToMatPoint() - lastRotation * (oldPoint.ToMatPoint());

    cv::Mat finalTransform = ImplantTools::JoinRigidTransform(lastRotation, lastTranslation) * mFemurTransformImplantToBone;

    return ImplantTools::getITKTransformFromCV(finalTransform);
}

itk::Rigid3DTransform<>::Pointer Balance::getNewImplantToBoneTibiaTransformAxial(const std::vector<PointTypeITK>& axialPointsCT)
{
    bool result;
    Plane axial = Plane::getBestPlane(axialPointsCT, result);
    axial.reverseByPoint(knee_.getAnkleCenter(), false);

    Point newNormal = axial.getNormalVector();
    Point oldNormal = PlaneTibia.getNormalVector();
    Point oldPoint = PlaneTibia.getPoint();

    cv::Mat lastRotation = ImplantTools::GetGeneralRotateTransformVectors(oldNormal, newNormal);
    cv::Mat lastTranslation = axial.getPoint().ToMatPoint() - lastRotation * (oldPoint.ToMatPoint());

    cv::Mat finalTransform = ImplantTools::JoinRigidTransform(lastRotation, lastTranslation) * mTibiaTransformImplantToBone;

    return ImplantTools::getITKTransformFromCV(finalTransform);
}

itk::Rigid3DTransform<>::Pointer Balance::getNewImplantToBoneFemurTransformCoronal(const std::vector<PointTypeITK>& coronalPointsCT)
{
    bool result;
    Plane coronal = Plane::getBestPlane(coronalPointsCT, result);
    Point posteriorPoint = knee_.getFemurKneeCenter() + 1000 * knee_.getFemurDirectVectorAP();
    coronal.reverseByPoint(posteriorPoint);

    Point newNormal = coronal.getNormalVector();
    Point oldNormal = PlaneA.getNormalVector();
    Point oldPoint = PlaneA.getPoint();

    cv::Mat lastRotation = ImplantTools::GetGeneralRotateTransformVectors(oldNormal, newNormal);
    cv::Mat lastTranslation = coronal.getPoint().ToMatPoint() - lastRotation * (oldPoint.ToMatPoint());

    cv::Mat finalTransform = ImplantTools::JoinRigidTransform(lastRotation, lastTranslation) * mFemurTransformImplantToBone;

    return ImplantTools::getITKTransformFromCV(finalTransform);
}
