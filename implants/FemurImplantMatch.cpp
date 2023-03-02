
#include "FemurImplantMatch.hpp"
#include <fstream>
#include "ImplantsException.hpp"
#include "ConvexHull.hpp"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkPlaneCollection.h"
#include "vtkClipClosedSurface.h"
#include "ImplantTools.hpp"

struct PointsBorder
{
    double latDistance = -1;
    double medDistance = -1;
    Point latNearPoint, medNearPoint;
    std::vector<Point> latPoints, medPoints;
};

inline cv::Mat Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer pTransform)
{
    itk::Matrix< double, 3, 3 > rotation = pTransform->GetMatrix();
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

inline cv::Mat Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer pTransform)
{
    itk::Vector< double, 3 > translate = pTransform->GetOffset();
    cv::Mat result(3, 1, CV_64FC1);

    result.at <double>(0, 0) = translate[0];
    result.at <double>(1, 0) = translate[1];
    result.at <double>(2, 0) = translate[2];
    return result;
}

inline cv::Mat GetRotateTransformAxisZ(Plane myPlane)
{
    //myPlane.normalizeNormalVector();
    Point normalXY(0, 0, 1);
    Point rotationAxis = normalXY.cross(myPlane.getNormalVector());
    rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
    double rotationAngle = Line::getAngleBetweenVectors(normalXY, myPlane.getNormalVector());
    cv::Mat rotation_1 = Line::getRotateMatrix(rotationAxis, -rotationAngle);
    cv::Mat rotation_2 = Line::getRotateMatrix(rotationAxis, rotationAngle);
    cv::Mat rotateVector_1 = rotation_1 * myPlane.getNormalVectorMat();
    cv::Mat rotateVector_2 = rotation_2 * myPlane.getNormalVectorMat();

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

    return rotate;
}

inline cv::Mat GetRotateTransformAxisX(Point sagitalVector, const cv::Mat& rotationZ)
{
    sagitalVector.normalice();
    cv::Mat finalSagitalMat = rotationZ * sagitalVector.ToMatPoint();
    Point finalSagitalVector = Point(finalSagitalMat);
    Point normalXY(1, 0, 0);
    Point rotationAxis = normalXY.cross(finalSagitalVector);
    rotationAxis = rotationAxis / sqrt(rotationAxis.dot(rotationAxis));
    double rotationAngle = Line::getAngleBetweenVectors(normalXY, finalSagitalVector);
    cv::Mat rotation_1 = Line::getRotateMatrix(rotationAxis, -rotationAngle);
    cv::Mat rotation_2 = Line::getRotateMatrix(rotationAxis, rotationAngle);
    cv::Mat rotateVector_1 = rotation_1 * finalSagitalMat;
    cv::Mat rotateVector_2 = rotation_2 * finalSagitalMat;

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

    return rotate;
}
inline cv::Mat getTransformToRobot(Plane currentPlane, const Plane& sagitalAnatomicPlane, const Point& P1, const Point& P2)
{
    Point sagitalP1 = sagitalAnatomicPlane.getProjectionPoint(P1);
    Point sagitalP2 = sagitalAnatomicPlane.getProjectionPoint(P2);
    sagitalP1 = currentPlane.getProjectionPoint(sagitalP1);
    sagitalP2 = currentPlane.getProjectionPoint(sagitalP2);
    Point sagitalVector = sagitalP2 - sagitalP1;

    //Se rota el plano hacia el eje Z.
    cv::Mat rotateZ = GetRotateTransformAxisZ(currentPlane);

    //Una vez el plano es perpendicular al eje Z, se gira hasta que su vector sagital sea paralelo al eje X.

    cv::Mat rotateX = GetRotateTransformAxisX(sagitalVector, rotateZ);
    cv::Mat finalMatrix = rotateX * rotateZ;
    return finalMatrix;
}


Plane FemurImplantMatch::finalTransformPlane(const Plane& plane, const itk::Rigid3DTransform<>::Pointer pTransform) const
{
    cv::Mat rotation = Rigid3DTransformToCVRotation(pTransform);
    cv::Mat translation = Rigid3DTransformToCVTranslation(pTransform);

    cv::Mat transformNormalVector = rotation * plane.getNormalVectorMat();
    cv::Mat transformPoint = (rotation * plane.getPointMat()) + translation;
    Plane transformPlane;
    transformPlane.init(Point(transformNormalVector), Point(transformPoint));
    return transformPlane;
}

FemurImplantMatch::FemurImplantMatch()
{
    isInit = false;
}
void FemurImplantMatch::init(const FemurImplant& implant, const Knee& knee, bool useKneeCenterAlignment)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_FEMUR_IMPLANT_MATCH;
    }
    this->implant = implant;
    this->knee = knee;
	this->useKneeCenterAlignment = useKneeCenterAlignment;
    getRotationMatrix();
    bool result = getTranslationMatrix();
    if (result == false)
    {
        throw ImplantExceptionCode::FAILED_TRANSFORMATION_MATCH_BY_CONDYLE;
    }

    result = getTranslationMatrixByCortex();
    if (result == false)
    {
        throw ImplantExceptionCode::FAILED_TRANSFORMATION_MATCH_BY_CORTEX;
    }
    isInit = true;
}

void FemurImplantMatch::getRotationMatrix()
{
    std::vector<cv::Point3d> implantVectors;
    std::vector<cv::Point3d> kneeVectors;

    Point implantTEA = implant.getDirectVectorTEA();
    Point implantFemurAxis = implant.getDirectVectorFemurAxis();
    Point implantAP = implant.getDirectVectorAP();
    implantVectors.push_back(implantTEA.ToCVPoint());
    implantVectors.push_back(implantFemurAxis.ToCVPoint());
    implantVectors.push_back(implantAP.ToCVPoint());

    Point kneeTEA = knee.getFemurVectorTEA();
    Point kneeFemurAxis = knee.getDirectVectorFemurAxis();
    Point kneeAP = knee.getFemurDirectVectorAP();

    kneeVectors.push_back(kneeTEA.ToCVPoint());
    kneeVectors.push_back(kneeFemurAxis.ToCVPoint());
    kneeVectors.push_back(kneeAP.ToCVPoint());

    cv::Mat implantMatrix = cv::Mat(implantVectors.size(), 3, CV_64F, implantVectors.data());
    cv::Mat kneeMatrix = cv::Mat(kneeVectors.size(), 3, CV_64F, kneeVectors.data());

    cv::Mat inverse = (implantMatrix.t()).inv();
    rotationMatrix = (kneeMatrix.t()) * inverse;
    //std::cout << "Rotation: " << rotationMatrix << std::endl;
}

bool FemurImplantMatch::getTranslationMatrix()
{
	Point tCenter;
	if (useKneeCenterAlignment = true)
	{
		tCenter = knee.getFemurKneeCenter();
	}
	else
	{
		tCenter = (knee.getLateralEpicondyle() + knee.getMedialEpicondyle()) / 2.;
	}

    cv::Mat kneeNormalVectorPlaneA = rotationMatrix * implant.getPlaneA().getNormalVectorMat();
    Plane kneePlaneA;
    kneePlaneA.init(Point(kneeNormalVectorPlaneA), knee.getMoveCondyle(implant.getImplantInfo()));

    cv::Mat kneeNormalVectorMidPlane = rotationMatrix * implant.getMidPlane().getNormalVectorMat();
    Plane kneeMidPlane;
    kneeMidPlane.init(Point(kneeNormalVectorMidPlane), tCenter);

    cv::Mat kneeNormalVectorPlaneC = rotationMatrix * implant.getPlaneC().getNormalVectorMat();
    Plane kneePlaneC;
    kneePlaneC.init(Point(kneeNormalVectorPlaneC), knee.getInferiorMoveFemurPoint(implant.getImplantInfo()));

    cv::Mat pSeudoExpKneePointA = rotationMatrix * implant.getPlaneA().getPointMat();
    cv::Mat pSeudoExpKneeMidPoint = rotationMatrix * implant.getMidPlane().getPointMat();
    cv::Mat pSeudoExpKneePointC = rotationMatrix * implant.getPlaneC().getPointMat();

    Point pSeudoKneePointA(pSeudoExpKneePointA);
    Point pSeudoKneeMidPoint(pSeudoExpKneeMidPoint);
    Point pSeudoKneePointC(pSeudoExpKneePointC);

    std::vector<cv::Point3d> normalVectors;
    std::vector<double> biasVector;
    double bias = 0.0;
    normalVectors.push_back(kneePlaneA.getNormalVector().ToCVPoint());
    normalVectors.push_back(kneeMidPlane.getNormalVector().ToCVPoint());
    normalVectors.push_back(kneePlaneC.getNormalVector().ToCVPoint());
    cv::Mat A(normalVectors.size(), 3, CV_64F, normalVectors.data());
    bias = -(kneePlaneA.getBias() + (pSeudoKneePointA.dot(kneePlaneA.getNormalVector())));
    biasVector.push_back(bias);
    bias = -(kneeMidPlane.getBias() + (pSeudoKneeMidPoint.dot(kneeMidPlane.getNormalVector())));
    biasVector.push_back(bias);
    bias = -(kneePlaneC.getBias() + (pSeudoKneePointC.dot(kneePlaneC.getNormalVector())));
    biasVector.push_back(bias);
    cv::Mat B(biasVector.size(), 1, CV_64F, biasVector.data());
    bool result = cv::solve(A, B, translationMatrix);
    //std::cout << "Translation: " << translationMatrix << std::endl;
    return result;
}

bool FemurImplantMatch::getTranslationMatrixByCortex()
{
	Point tCenter;
	if (useKneeCenterAlignment = true)
	{
		tCenter = knee.getFemurKneeCenter();
	}
	else
	{
		tCenter = (knee.getLateralEpicondyle() + knee.getMedialEpicondyle()) / 2.;
	}

    cv::Mat kneeNormalVectorPlaneE = rotationMatrix * implant.getPlaneE().getNormalVectorMat();
    Plane kneePlaneE;
    kneePlaneE.init(Point(kneeNormalVectorPlaneE), knee.getAnteriorCortex());

    cv::Mat kneeNormalVectorMidPlane = rotationMatrix * implant.getMidPlane().getNormalVectorMat();
    Plane kneeMidPlane;
    kneeMidPlane.init(Point(kneeNormalVectorMidPlane), tCenter);

    cv::Mat kneeNormalVectorPlaneC = rotationMatrix * implant.getPlaneC().getNormalVectorMat();
    Plane kneePlaneC;
    kneePlaneC.init(Point(kneeNormalVectorPlaneC), knee.getInferiorMoveFemurPoint(implant.getImplantInfo()));

    cv::Mat pSeudoExpKneePointE = rotationMatrix * implant.getPlaneE().getPointMat();
    cv::Mat pSeudoExpKneeMidPoint = rotationMatrix * implant.getMidPlane().getPointMat();
    cv::Mat pSeudoExpKneePointC = rotationMatrix * implant.getPlaneC().getPointMat();

    Point pSeudoKneePointE(pSeudoExpKneePointE);
    Point pSeudoKneeMidPoint(pSeudoExpKneeMidPoint);
    Point pSeudoKneePointC(pSeudoExpKneePointC);

    std::vector<cv::Point3d> normalVectors;
    std::vector<double> biasVector;
    double bias = 0.0;
    normalVectors.push_back(kneePlaneE.getNormalVector().ToCVPoint());
    normalVectors.push_back(kneeMidPlane.getNormalVector().ToCVPoint());
    normalVectors.push_back(kneePlaneC.getNormalVector().ToCVPoint());
    cv::Mat A(normalVectors.size(), 3, CV_64F, normalVectors.data());
    bias = -(kneePlaneE.getBias() + (pSeudoKneePointE.dot(kneePlaneE.getNormalVector())));
    biasVector.push_back(bias);
    bias = -(kneeMidPlane.getBias() + (pSeudoKneeMidPoint.dot(kneeMidPlane.getNormalVector())));
    biasVector.push_back(bias);
    bias = -(kneePlaneC.getBias() + (pSeudoKneePointC.dot(kneePlaneC.getNormalVector())));
    biasVector.push_back(bias);
    cv::Mat B(biasVector.size(), 1, CV_64F, biasVector.data());
    bool result = cv::solve(A, B, translationMatrixByCortex);
    //std::cout << "Translation 2: " << translationMatrixByCortex << std::endl;
    return result;
}

vtkSmartPointer<vtkPolyData> FemurImplantMatch::GetCuttingFemur(bool translateByCondyle)
{
    Point anteriorPoint, posteriorPoint;
    double resizeVector = 1000000.0;

    anteriorPoint = knee.getFemurKneeCenter() + resizeVector * knee.getFemurDirectVectorAP();
    posteriorPoint = knee.getFemurKneeCenter() - resizeVector * knee.getFemurDirectVectorAP();

    Plane PlaneA = transformPlane(implant.getPlaneA(), translateByCondyle);
    PlaneA.reverseByPoint(anteriorPoint);

    Plane PlaneB = transformPlane(implant.getPlaneB(), translateByCondyle);
    PlaneB.reverseByPoint(anteriorPoint);

    Plane PlaneC = transformPlane(implant.getPlaneC(), translateByCondyle);
    PlaneC.reverseByPoint(knee.getHipCenter());

    Plane PlaneD = transformPlane(implant.getPlaneD(), translateByCondyle);
    PlaneD.reverseByPoint(posteriorPoint);

    Plane PlaneE = transformPlane(implant.getPlaneE(), translateByCondyle);
    PlaneE.reverseByPoint(posteriorPoint);

    cv::Point3d planeNormal, planePoint;

    vtkNew<vtkPlane> vtkPlaneA, vtkPlaneB, vtkPlaneC, vtkPlaneD, vtkPlaneE;

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

    planeNormal = PlaneE.getNormalVector();
    planePoint = PlaneE.getPoint();
    vtkPlaneE->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
    vtkPlaneE->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);

    vtkNew<vtkPlaneCollection> femurPlanes;
    femurPlanes->AddItem(vtkPlaneA);
    femurPlanes->AddItem(vtkPlaneB);
    femurPlanes->AddItem(vtkPlaneC);
    femurPlanes->AddItem(vtkPlaneD);
    femurPlanes->AddItem(vtkPlaneE);

    vtkNew<vtkClipClosedSurface> femurClipper;
    femurClipper->SetInputData(knee.GetFemurPoly());
    femurClipper->SetClippingPlanes(femurPlanes);
    femurClipper->Update();

    return femurClipper->GetOutput();
}

Plane FemurImplantMatch::transformPlane(const Plane& plane, bool translateByCondyle) const
{
    cv::Mat transformNormalVector = rotationMatrix * plane.getNormalVectorMat();
    cv::Mat transformPoint;
    if (translateByCondyle == true)
        transformPoint = (rotationMatrix * plane.getPointMat()) + translationMatrix;
    else
        transformPoint = (rotationMatrix * plane.getPointMat()) + translationMatrixByCortex;
    Plane transformPlane;
    transformPlane.init(Point(transformNormalVector), Point(transformPoint));
    return transformPlane;
}

Plane FemurImplantMatch::GetPlane(PlaneID id, bool translateByCondyle) const
{
    if (id == kPlaneA)
    {
        return transformPlane(implant.getPlaneA(), translateByCondyle);
    }
    else if (id == kPlaneB)
    {
        return transformPlane(implant.getPlaneB(), translateByCondyle);
    }
    else if (id == kPlaneC)
    {
        return transformPlane(implant.getPlaneC(), translateByCondyle);
    }
    else if (id == kPlaneD)
    {
        return transformPlane(implant.getPlaneD(), translateByCondyle);
    }
    else if (id == kPlaneE)
    {
        return transformPlane(implant.getPlaneE(), translateByCondyle);
    }
    else if (id == kPlaneMid)
    {
        return transformPlane(implant.getMidPlane(), translateByCondyle);
    }
    else
    {
        return Plane();
    }
}

Point FemurImplantMatch::TransformImplantPointToBone(const Point& pPoint, bool translateByCondyle) const
{
    cv::Mat mat(3, 1, CV_64FC1);
    mat.at <double>(0, 0) = pPoint.x;
    mat.at <double>(1, 0) = pPoint.y;
    mat.at <double>(2, 0) = pPoint.z;

    cv::Mat transformPoint;
    if (translateByCondyle == true)
        transformPoint = (rotationMatrix * mat) + translationMatrix;
    else
        transformPoint = (rotationMatrix * mat) + translationMatrixByCortex;
    return Point(transformPoint);
}

Knee FemurImplantMatch::getKnee() const
{
    return knee;
}

vtkSmartPointer<vtkPolyData> FemurImplantMatch::getContour(const vtkSmartPointer<vtkPolyData> poly, const Point& pNormal, const Point& pPoint) const
{
    vtkNew<vtkPlane> cutPlane;
    cutPlane->SetOrigin(pPoint.x, pPoint.y, pPoint.z);
    cutPlane->SetNormal(pNormal.x, pNormal.y, pNormal.z);

    vtkNew<vtkCutter> cutter;
    cutter->SetInputData(poly);
    cutter->SetCutFunction(cutPlane);
    cutter->Update();

    auto contour = cutter->GetOutput();

    return contour;
}

Point FemurImplantMatch::getNearPointUnderCortex(const Plane& myPlane, std::vector<Point>& points, double distance) const
{
    Point cortex = knee.getAnteriorCortex();
    Point hip = knee.getHipCenter();
    Point femurKnee = knee.getFemurKneeCenter();
    Point directVector = hip - femurKnee;
    Plane Transverse;
    Transverse.init(directVector, cortex);
    int transverseSign;
    if (Transverse.eval(hip) > 0)
        transverseSign = -1.0;
    else
        transverseSign = 1.0;

    std::vector<Point> femurPoints = knee.getFemurPoints();

    auto it1 = femurPoints.begin();
    auto it2 = femurPoints.end();

    Point midPoint, tempPoint;
    points.clear();

    if (distance == 0)
    {
        vtkSmartPointer<vtkPolyData> contour = getContour(knee.GetFemurPoly(), myPlane.getNormalVector(), myPlane.getPoint());
        vtkSmartPointer<vtkPoints> pointsCut = contour->GetPoints();
        int tSize = pointsCut->GetNumberOfPoints();
        for (int i = 0; i < tSize; i++)
        {
            double pnt[3];
            pointsCut->GetPoint(i, pnt);
            Point myPoint(pnt[0], pnt[1], pnt[2]);

            if (transverseSign * Transverse.eval(myPoint) > 0)
            {
                midPoint = midPoint + myPoint;
                points.push_back(myPoint);
            }
        }

        if (points.size() > 0)
        {
            midPoint = midPoint / double(points.size());
        }
    }
    else
    {
        for (; it1 != it2; ++it1)
        {
            if (transverseSign * Transverse.eval(*it1) > 0 && myPlane.isPointNearToPlane(*it1, distance) == true)
            {
                tempPoint = myPlane.getProjectionPoint(*it1);
                midPoint = midPoint + tempPoint;
                points.push_back(tempPoint);
            }
        }

        if (points.size() > 0)
        {
            midPoint = midPoint / double(points.size());
        }
    }

    return midPoint;
}

std::vector<Point> FemurImplantMatch::GetPointsNearPlane(PlaneID id, bool translateByCondyle, double distance) const
{
    std::vector<Point> temp;

    if (id == kPlaneA)
    {
        Plane plane = transformPlane(implant.getPlaneA(), translateByCondyle);
        getNearPointUnderCortex(plane, temp, distance);
    }
    else if (id == kPlaneB)
    {
        Plane plane = transformPlane(implant.getPlaneB(), translateByCondyle);
        getNearPointUnderCortex(plane, temp, distance);
    }
    else if (id == kPlaneC)
    {
        Plane plane = transformPlane(implant.getPlaneC(), translateByCondyle);
        getNearPointUnderCortex(plane, temp, distance);
    }
    else if (id == kPlaneD)
    {
        Plane plane = transformPlane(implant.getPlaneD(), translateByCondyle);
        getNearPointUnderCortex(plane, temp, distance);
    }
    else if (id == kPlaneE)
    {
        Plane plane = transformPlane(implant.getPlaneE(), translateByCondyle);
        getNearPointUnderCortex(plane, temp, distance);
    }
    else
    {
        return temp;
    }

    return temp;
}

itk::Matrix< double, 3, 3 > FemurImplantMatch::GetRotationMatrix() const
{
    itk::Matrix< double, 3, 3 > rotation;

    rotation[0][0] = rotationMatrix.at <double>(0, 0);
    rotation[0][1] = rotationMatrix.at <double>(0, 1);
    rotation[0][2] = rotationMatrix.at <double>(0, 2);

    rotation[1][0] = rotationMatrix.at <double>(1, 0);
    rotation[1][1] = rotationMatrix.at <double>(1, 1);
    rotation[1][2] = rotationMatrix.at <double>(1, 2);

    rotation[2][0] = rotationMatrix.at <double>(2, 0);
    rotation[2][1] = rotationMatrix.at <double>(2, 1);
    rotation[2][2] = rotationMatrix.at <double>(2, 2);
    return rotation;
}

itk::Vector< double, 3 > FemurImplantMatch::GetTranslationMatrix() const
{
    itk::Vector< double, 3 > translation;
    translation[0] = translationMatrix.at <double>(0, 0);
    translation[1] = translationMatrix.at <double>(1, 0);
    translation[2] = translationMatrix.at <double>(2, 0);

    return translation;
}

itk::Vector< double, 3 > FemurImplantMatch::GetTranslationMatrixByCortex() const
{
    itk::Vector< double, 3 > translation;
    translation[0] = translationMatrixByCortex.at <double>(0, 0);
    translation[1] = translationMatrixByCortex.at <double>(1, 0);
    translation[2] = translationMatrixByCortex.at <double>(2, 0);

    return translation;
}

void savePoints(std::vector<Point> points, cv::Mat myRotation, std::string name)
{
    std::ofstream outfilePoint(name + ".txt");
    for (int i = 0; i < points.size(); i++)
    {
        Point tempPointPrint = points[i];
        cv::Mat tempPointRotMat = myRotation * (tempPointPrint).ToMatPoint();
        Point printPoint = Point(tempPointRotMat);
        outfilePoint << printPoint.x << " " << printPoint.y << "\n";
    }
    outfilePoint.close();
}

std::vector<PointTypeITK> FemurImplantMatch::GetHullPoints(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut, PlaneID id, double distanceSide, double distanceTop, double angleLateral, double angleMedial, int amount, bool posteriorLongCurve) const
{
    std::vector<Point> projectedPoints;
    Point midPointPlane;
    std::vector<PointTypeITK> hull;
    double resizeVector = 1000000.0;

    Point centerP1, centerP2, topPoint, downPoint, lateralPoint, medialPoint, lateralSide, medialSide;
    Plane midPlane = finalTransformPlane(implant.getMidPlane(), pTransformIn);
    Point fromPostToAntVector = knee.getFemurDirectVectorAP();

    lateralPoint = knee.getLateralEpicondyle();

    medialPoint = knee.getMedialEpicondylePerp();

    Point fromMedToLatVector = lateralPoint - medialPoint;
    fromMedToLatVector = fromMedToLatVector / sqrt(fromMedToLatVector.dot(fromMedToLatVector));

    Point anterior = knee.getFemurKneeCenter() + (resizeVector * fromPostToAntVector);

    Point posterior = knee.getFemurKneeCenter() - (resizeVector * fromPostToAntVector);

    lateralPoint = knee.getFemurKneeCenter() + (resizeVector * fromMedToLatVector);

    medialPoint = knee.getFemurKneeCenter() - (resizeVector * fromMedToLatVector);
    Plane currentPlane, sagitalAnatomicPlane;
    sagitalAnatomicPlane.init(knee.getFemurVectorTEA(), knee.getFemurKneeCenter());
    std::string planeName;
    cv::Mat myRotation;
    Point myNormal, myNormalTemp;
    if (id == kPlaneA)
    {
        planeName = "points_A";
        currentPlane = finalTransformPlane(implant.getPlaneA(), pTransformIn);
        projectedPoints.clear();
        midPointPlane = getNearPointUnderCortex(currentPlane, projectedPoints);

        myNormalTemp = knee.getFemurKneeCenter() - resizeVector * knee.getFemurDirectVectorAP();
        myNormal = currentPlane.getProjectionPoint(myNormalTemp) - myNormalTemp;

        currentPlane.fixNormalVector(myNormal);

        if (projectedPoints.size() > 15)
        {
            centerP1 = midPointPlane;
            centerP2 = currentPlane.getProjectionPoint(knee.getHipCenter());

            myRotation = getTransformToRobot(currentPlane, sagitalAnatomicPlane, centerP1, centerP2);

            centerP1 = midPlane.getProjectionPoint(centerP1);
            centerP2 = midPlane.getProjectionPoint(centerP2);

            topPoint = currentPlane.getProjectionPoint(knee.getHipCenter());
            downPoint = currentPlane.getProjectionPoint(knee.getAnkleCenter());
            lateralSide = currentPlane.getProjectionPoint(lateralPoint);
            medialSide = currentPlane.getProjectionPoint(medialPoint);

        }
        else
        {
            return hull;
        }
    }
    else if (id == kPlaneB)
    {
        planeName = "points_B";
        currentPlane = finalTransformPlane(implant.getPlaneB(), pTransformIn);
        projectedPoints.clear();
        midPointPlane = getNearPointUnderCortex(currentPlane, projectedPoints);

        myNormalTemp = knee.getFemurKneeCenter() - resizeVector * knee.getFemurDirectVectorAP();
        myNormal = myNormalTemp - currentPlane.getProjectionPoint(myNormalTemp);

        currentPlane.fixNormalVector(myNormal);

        if (projectedPoints.size() > 15)
        {
            centerP1 = midPointPlane;
            centerP2 = currentPlane.getProjectionPoint(knee.getHipCenter());

            myRotation = getTransformToRobot(currentPlane, sagitalAnatomicPlane, centerP1, centerP2);

            centerP1 = midPlane.getProjectionPoint(centerP1);
            centerP2 = midPlane.getProjectionPoint(centerP2);

            topPoint = currentPlane.getProjectionPoint(posterior);
            downPoint = currentPlane.getProjectionPoint(anterior);
            lateralSide = currentPlane.getProjectionPoint(lateralPoint);
            medialSide = currentPlane.getProjectionPoint(medialPoint);


        }
        else
        {
            return hull;
        }
    }
    else if (id == kPlaneC)
    {
        planeName = "points_C";
        currentPlane = finalTransformPlane(implant.getPlaneC(), pTransformIn);
        projectedPoints.clear();
        midPointPlane = getNearPointUnderCortex(currentPlane, projectedPoints);

        myNormalTemp = knee.getAnkleCenter();
        myNormal = myNormalTemp - currentPlane.getProjectionPoint(myNormalTemp);

        currentPlane.fixNormalVector(myNormal);

        if (projectedPoints.size() > 15)
        {
            centerP1 = midPointPlane;
            centerP2 = knee.getFemurKneeCenter() - resizeVector * knee.getFemurDirectVectorAP();
            centerP2 = currentPlane.getProjectionPoint(centerP2);

            myRotation = getTransformToRobot(currentPlane, sagitalAnatomicPlane, centerP1, centerP2);

            centerP1 = midPlane.getProjectionPoint(centerP1);
            centerP2 = midPlane.getProjectionPoint(centerP2);

            topPoint = currentPlane.getProjectionPoint(posterior);
            downPoint = currentPlane.getProjectionPoint(anterior);
            lateralSide = currentPlane.getProjectionPoint(lateralPoint);
            medialSide = currentPlane.getProjectionPoint(medialPoint);
        }
        else
        {
            return hull;
        }
    }
    else if (id == kPlaneD)
    {
        planeName = "points_D";
        currentPlane = finalTransformPlane(implant.getPlaneD(), pTransformIn);
        projectedPoints.clear();
        midPointPlane = getNearPointUnderCortex(currentPlane, projectedPoints);

        myNormalTemp = knee.getFemurKneeCenter() + resizeVector * knee.getFemurDirectVectorAP();
        myNormal = myNormalTemp - currentPlane.getProjectionPoint(myNormalTemp);

        currentPlane.fixNormalVector(myNormal);

        if (projectedPoints.size() > 15)
        {
            centerP1 = midPointPlane;
            centerP2 = currentPlane.getProjectionPoint(knee.getHipCenter());

            myRotation = getTransformToRobot(currentPlane, sagitalAnatomicPlane, centerP1, centerP2);

            centerP1 = midPlane.getProjectionPoint(centerP1);
            centerP2 = midPlane.getProjectionPoint(centerP2);

            topPoint = currentPlane.getProjectionPoint(anterior);
            downPoint = currentPlane.getProjectionPoint(posterior);
            lateralSide = currentPlane.getProjectionPoint(lateralPoint);
            medialSide = currentPlane.getProjectionPoint(medialPoint);

        }
        else
        {
            return hull;
        }
    }
    else if (id == kPlaneE)
    {
        planeName = "points_E";
        currentPlane = finalTransformPlane(implant.getPlaneE(), pTransformIn);
        projectedPoints.clear();
        midPointPlane = getNearPointUnderCortex(currentPlane, projectedPoints);

        myNormalTemp = knee.getFemurKneeCenter() + resizeVector * knee.getFemurDirectVectorAP();
        myNormal = myNormalTemp - currentPlane.getProjectionPoint(myNormalTemp);

        currentPlane.fixNormalVector(myNormal);

        if (projectedPoints.size() > 15)
        {
            centerP1 = midPointPlane;
            centerP2 = currentPlane.getProjectionPoint(knee.getHipCenter());

            myRotation = getTransformToRobot(currentPlane, sagitalAnatomicPlane, centerP1, centerP2);

            centerP1 = midPlane.getProjectionPoint(centerP1);
            centerP2 = midPlane.getProjectionPoint(centerP2);

            topPoint = currentPlane.getProjectionPoint(knee.getHipCenter());
            downPoint = currentPlane.getProjectionPoint(knee.getAnkleCenter());
            lateralSide = currentPlane.getProjectionPoint(lateralPoint);
            medialSide = currentPlane.getProjectionPoint(medialPoint);

        }
        else
        {
            return hull;
        }
    }
    else
    {
        return hull;
    }
    std::vector<Point> vertices, cutPoints;

    if (angleLateral > 45)
    {
        angleLateral = 45;
    }

    if (angleLateral < 0)
    {
        angleLateral = 0;
    }

    if (angleMedial > 45)
    {
        angleMedial = 45;
    }

    if (angleMedial < 0)
    {
        angleMedial = 0;
    }

    double angleLatRad = ((90.0 - angleLateral) * PI) / 180.0;
    double angleMedRad = ((90.0 - angleMedial) * PI) / 180.0;

    if (id == kPlaneC || id == kPlaneD || id == kPlaneE)
    {
        getVerticesCDE(projectedPoints, downPoint, lateralSide, medialSide, topPoint, centerP1, centerP2, distanceSide, distanceTop, angleLatRad, angleMedRad, vertices);
    }
    else if (id == kPlaneA)
    {
        getVerticesA(projectedPoints, downPoint, lateralSide, medialSide, topPoint, midPlane, currentPlane, myRotation, vertices, distanceSide, amount, posteriorLongCurve);
    }
    else
    {
        getVerticesB(projectedPoints, downPoint, lateralSide, medialSide, topPoint, midPlane, currentPlane, myRotation, vertices, distanceSide, distanceTop, amount);
    }

    Point initExtreme, endExtreme;

    if (vertices.size() > 1)
    {
        initExtreme = vertices[0];
        endExtreme = vertices[vertices.size() - 1];

        cv::Mat initExtremeMat, endExtremeMat;

        initExtremeMat = myRotation * initExtreme.ToMatPoint();
        endExtremeMat = myRotation * endExtreme.ToMatPoint();

        initExtreme = Point(initExtremeMat);
        endExtreme = Point(endExtremeMat);

        if (endExtreme.y < initExtreme.y)
        {
            std::reverse(vertices.begin(), vertices.end());
        }
    }

    hull = increaseVectorToAmount(vertices, amount);

    itk::Vector< double, 3 > translate;
    if (projectedPoints.size() > 0)
    {
        midPointPlane = sagitalAnatomicPlane.getProjectionPoint(midPointPlane);
        midPointPlane = currentPlane.getProjectionPoint(midPointPlane);
        Point tTemp = myRotation * (midPointPlane.ToMatPoint());
        translate[0] = -tTemp.x;
        translate[1] = -tTemp.y;
        translate[2] = -tTemp.z;
    }
    else
    {
        translate[0] = 0;
        translate[1] = 0;
        translate[2] = 0;
    }

    itk::Matrix< double, 3, 3 > rotationITK;
    rotationITK[0][0] = myRotation.at<double>(0, 0);
    rotationITK[0][1] = myRotation.at<double>(0, 1);
    rotationITK[0][2] = myRotation.at<double>(0, 2);

    rotationITK[1][0] = myRotation.at<double>(1, 0);
    rotationITK[1][1] = myRotation.at<double>(1, 1);
    rotationITK[1][2] = myRotation.at<double>(1, 2);

    rotationITK[2][0] = myRotation.at<double>(2, 0);
    rotationITK[2][1] = myRotation.at<double>(2, 1);
    rotationITK[2][2] = myRotation.at<double>(2, 2);

    pTransformOut->SetMatrix(rotationITK);
    pTransformOut->SetOffset(translate);


    //////////////////////////////////////////////////////////////////////////

    //std::ofstream outfilePoint(planeName + ".txt");
    //std::ofstream outfileBorder(planeName + "_border.txt");
    //for (int i = 0; i < hull.size(); i++)
    //{
    //    Point tempPointPrint(hull[i][0], hull[i][1], hull[i][2]);
    //    //std::cout << tempPointPrint << std::endl;
    //    cv::Mat tempPointRotMat = myRotation * (tempPointPrint).ToMatPoint();
    //    Point printPoint = Point(tempPointRotMat);
    //    outfileBorder << printPoint.x << " " << printPoint.y << "\n";
    //}

    //for (int i = 0; i < projectedPoints.size(); i++)
    //{
    //    cv::Mat tempPointRotMat = myRotation * (projectedPoints[i]).ToMatPoint();
    //    Point printPoint = Point(tempPointRotMat);
    //    outfilePoint << printPoint.x << " " << printPoint.y << "\n";
    //}
    //outfilePoint.close();
    //outfileBorder.close();

    ////////////////////////////////////////////////////////////////////////////////////////////////


    //std::ofstream outfile(planeName);
    //double a, b, c;
    //std::vector<Point> printPoints;
    //if (vertices.size() > 0)
    //{
    //    printPoints = vertices;
    //    for (int i = 0; i < printPoints.size() - 1; i++)
    //    {
    //        //Point example1 = printPoints[i];
    //        Point P1 = printPoints[i];
    //        Point P2 = printPoints[i + 1];

    //        for (double j = 0; j < 100; j++)
    //        {
    //            Point example1 = (P1 + (j / 100)*(P2 - P1));
    //            //cv::Mat ppp = myRotation * example1.ToMatPoint();
    //            //example1 = Point(ppp);
    //            outfile << example1.x << " " << example1.y << " " << " " << example1.z << "\n";
    //        }
    //        //outfile << example1.x << " " << example1.y << " " << " " << example1.z << "\n";
    //    }
    //    outfile.close();
    //}
    return hull;
}

void FemurImplantMatch::getVerticesA(const std::vector<Point>& points, const Point& downPoint, const Point& lateralPoint, const Point& medialPoint, const Point& topPoint, const Plane& midPlane, const Plane& currentPlane, const cv::Mat& pRotation, std::vector<Point>& vertices, double distance, int amount, bool longCurve) const
{
    float maxDist = 15;
    if (abs(distance) > maxDist)
    {
        distance = (maxDist / float(distance)) * float(distance);
    }

    Plane myMidPlane = midPlane;
    myMidPlane.normalizeNormalVector();
    if (myMidPlane.eval(lateralPoint) < 0)
    {
        myMidPlane.reverse();
    }

    Point centerP1 = points[0];
    Point centerP2 = currentPlane.getProjectionPoint(knee.getHipCenter());

    Point vectorAP = myMidPlane.getProjectionPoint(centerP1) - myMidPlane.getProjectionPoint(centerP2);

    Line topLateral(midPlane.getNormalVector(), topPoint);
    Line topMedial(midPlane.getNormalVector(), topPoint);
    Line downLateral(midPlane.getNormalVector(), downPoint);
    Line downMedial(midPlane.getNormalVector(), downPoint);
    Line medialLine(vectorAP, medialPoint);
    Line lateralLine(vectorAP, lateralPoint);

    double latTopDistance = -1.0;
    double medTopDistance = -1.0;
    double latDownDistance = -1.0;
    double medDownDistance = -1.0;
    double medialDistance = -1.0;
    double lateralDistance = -1.0;
    double latTopTemp, medTopTemp, latDownTemp, medDownTemp, lateralTemp, medialTemp;
    Point latTopPoint, medTopPoint, latDownPoint, medDownPoint, latPoint, medPoint;
    std::vector<Point>::const_iterator it1, it2;
    it1 = points.begin();
    it2 = points.end();
    for (; it1 != it2; ++it1)
    {
        if (myMidPlane.eval(*it1) > 0)
        {
            latTopTemp = topLateral.getDistanceFromPoint(*it1);
            latDownTemp = downLateral.getDistanceFromPoint(*it1);
            lateralTemp = lateralLine.getDistanceFromPoint(*it1);

            if (latTopTemp > latTopDistance)
            {
                latTopDistance = latTopTemp;
                latDownPoint = *it1;
            }

            if (latDownTemp > latDownDistance)
            {
                latDownDistance = latDownTemp;
                latTopPoint = *it1;
            }

            if (lateralTemp > lateralDistance)
            {
                lateralDistance = lateralTemp;
                medPoint = *it1;
            }
        }
        else
        {
            medTopTemp = topMedial.getDistanceFromPoint(*it1);
            medDownTemp = downMedial.getDistanceFromPoint(*it1);
            medialTemp = medialLine.getDistanceFromPoint(*it1);

            if (medTopTemp > medTopDistance)
            {
                medTopDistance = medTopTemp;
                medDownPoint = *it1;
            }

            if (medDownTemp > medDownDistance)
            {
                medDownDistance = medDownTemp;
                medTopPoint = *it1;
            }

            if (medialTemp > medialDistance)
            {
                medialDistance = medialTemp;
                latPoint = *it1;
            }
        }
    }
    std::vector<Point> lateralOut, medialOut;
    Point centerPoint;
    separateMedialAndLateralPoints(points, midPlane, latPoint, medPoint, lateralOut, medialOut, centerPoint);

    if (lateralOut.size() == 0 || medialOut.size() == 0)
    {
        throw ImplantExceptionCode::CAN_NOT_DEFINE_DIFFERENT_BORDER_FOR_LATERAL_AND_MEDIAL_SIDE;
    }

    Plane sagitalPlane;
    sagitalPlane.init(midPlane.getNormalVector(), centerPoint);
    std::vector<Point> concaveLat, concaveMed, mainLatPoints, mainMedPoints;
    ConvexHull latHull(lateralOut, pRotation);
    ConvexHull medHull(medialOut, pRotation);
    concaveLat = latHull.GetConvexHull();
    concaveMed = medHull.GetConvexHull();

    if (concaveLat.size() == 0)
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_CONVEX_HULL_ON_LATERAL_SIDE;
    }

    if (concaveMed.size() == 0)
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_CONVEX_HULL_ON_MEDIAL_SIDE;
    }

    bool latChange, medChange;
    std::pair<int, int> posLat, posMed;

    latChange = getConcaveMainPoints(concaveLat, currentPlane, sagitalPlane, downPoint, 0, 1, mainLatPoints, posLat, 0.5);
    medChange = getConcaveMainPoints(concaveMed, currentPlane, sagitalPlane, downPoint, 0, 1, mainMedPoints, posMed, 0.5);

    Point refDown;
    Plane extremeDownPlane;
    extremeDownPlane.init(vectorAP, downPoint);
    if (extremeDownPlane.getDistanceFromPoint(medDownPoint) > extremeDownPlane.getDistanceFromPoint(latDownPoint))
    {
        refDown = longCurve ? medDownPoint : (medDownPoint + medTopPoint) / 2.;
    }
    else
    {
        refDown = longCurve ? latDownPoint : (latDownPoint + latTopPoint) / 2.;
    }

    extremeDownPlane.movePlane(refDown);
    Point centerDownPoint = extremeDownPlane.getProjectionPoint(centerPoint);
    bool isLatBegin = true;

    if (mainLatPoints.size() > 1)
    {
        //mainLatPoints.insert(mainLatPoints.begin(), extremeDownPlane.getProjectionPoint(mainLatPoints[0]));

        if (latChange == true)
        {
            std::reverse(mainLatPoints.begin(), mainLatPoints.end());
            isLatBegin = false;

            posLat.first = mainLatPoints.size() - posLat.first - 1;
            posLat.second = mainLatPoints.size() - posLat.second - 1;
        }
    }
    else
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_MAIN_POINTS_IN_CONVEX_HULL_ON_LATERAL_SIDE;
    }

    if (mainMedPoints.size() > 1)
    {
        //mainMedPoints.insert(mainMedPoints.begin(), extremeDownPlane.getProjectionPoint(mainMedPoints[0]));

        if (medChange == true)
        {
            std::reverse(mainMedPoints.begin(), mainMedPoints.end());

            posMed.first = mainMedPoints.size() - posMed.first - 1;
            posMed.second = mainMedPoints.size() - posMed.second - 1;
        }
    }
    else
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_MAIN_POINTS_IN_CONVEX_HULL_ON_MEDIAL_SIDE;
    }

    std::vector<Point> allVertices;

    if (isLatBegin == true)
    {
        for (int i = 0; i < mainLatPoints.size(); i++)
        {
            allVertices.push_back(mainLatPoints[i]);
        }

        allVertices.push_back(centerDownPoint);

        posMed.first = allVertices.size() + posMed.first;
        posMed.second = allVertices.size() + posMed.second;

        for (int i = 0; i < mainMedPoints.size(); i++)
        {
            allVertices.push_back(mainMedPoints[i]);
        }
    }
    else
    {
        for (int i = 0; i < mainMedPoints.size(); i++)
        {
            allVertices.push_back(mainMedPoints[i]);
        }

        allVertices.push_back(centerDownPoint);

        posLat.first = allVertices.size() + posLat.first;
        posLat.second = allVertices.size() + posLat.second;

        for (int i = 0; i < mainLatPoints.size(); i++)
        {
            allVertices.push_back(mainLatPoints[i]);
        }
    }

    int posInitOut, posEndOut, posInitIn, posEndIn;

    if (posLat.first < posMed.first)
    {
        posInitIn = posLat.first;
        posEndIn = posMed.first;
    }
    else
    {
        posInitIn = posMed.first;
        posEndIn = posLat.first;
    }

    if (posLat.second < posMed.second)
    {
        posInitOut = posLat.second;
        posEndOut = posMed.second;
    }
    else
    {
        posInitOut = posMed.second;
        posEndOut = posLat.second;
    }

    Point computePoint;

    for (int i = 0; i < allVertices.size() - 1; i++)
    {
        if (i >= posInitOut && i <= posEndOut)
        {
            vertices.push_back(allVertices[i]);
            continue;
        }

        if (distance != 0)
        {
            computePoint = movePointAtNormal(allVertices[i], allVertices[i + 1], pRotation, distance);
            vertices.push_back(computePoint);

            if (i == allVertices.size() - 2)
            {
                computePoint = movePointAtNormal(allVertices[i], allVertices[i + 1], pRotation, distance, true);
                vertices.push_back(computePoint);
            }
        }
        else
        {
            vertices.push_back(allVertices[i]);
        }
    }

    vertices = ConvexHull::interpolateSpline(vertices, amount);
}

void FemurImplantMatch::getVerticesB(const std::vector<Point>& points, const Point& downPoint, const Point& lateralPoint, const Point& medialPoint, const Point& topPoint, const Plane& midPlane, const Plane& currentPlane, const cv::Mat& pRotation, std::vector<Point>& vertices, double distanceSide, double distanceTop, int amount) const
{
    float maxDist = 15;

    if (abs(distanceSide) > maxDist)
    {
        distanceSide = (maxDist / float(distanceSide)) * float(distanceSide);
    }

    if (abs(distanceTop) > maxDist)
    {
        distanceTop = (maxDist / float(distanceTop)) * float(distanceTop);
    }

    Plane myMidPlane = midPlane;
    myMidPlane.normalizeNormalVector();
    if (myMidPlane.eval(lateralPoint) < 0)
    {
        myMidPlane.reverse();
    }

    Point centerP1 = points[0];
    Point centerP2 = currentPlane.getProjectionPoint(knee.getHipCenter());

    Point vectorAP = myMidPlane.getProjectionPoint(centerP1) - myMidPlane.getProjectionPoint(centerP2);

    Line downLateral(midPlane.getNormalVector(), downPoint);
    Line downMedial(midPlane.getNormalVector(), downPoint);
    Line medialLine(vectorAP, medialPoint);
    Line lateralLine(vectorAP, lateralPoint);

    double latTopDistance = -1.0;
    double medTopDistance = -1.0;
    double latDownDistance = -1.0;
    double medDownDistance = -1.0;
    double medialDistance = -1.0;
    double lateralDistance = -1.0;
    double latDownTemp, medDownTemp, lateralTemp, medialTemp;
    Point latTopPoint, medTopPoint, latPoint, medPoint, centerLatIn, centerMedIn, latDownPoint, medDownPoint;
    std::vector<Point>::const_iterator it1, it2;
    int contLat = 0, contMed = 0;
    it1 = points.begin();
    it2 = points.end();
    for (; it1 != it2; ++it1)
    {
        if (myMidPlane.eval(*it1) > 0)
        {
            latDownTemp = downLateral.getDistanceFromPoint(*it1);
            lateralTemp = lateralLine.getDistanceFromPoint(*it1);

            if (latDownTemp > latTopDistance)
            {
                latTopDistance = latDownTemp;
                latTopPoint = *it1;
            }

            if (lateralTemp > lateralDistance)
            {
                lateralDistance = lateralTemp;
                latPoint = *it1;
            }

            if (latDownTemp < latDownDistance || latDownDistance < 0)
            {
                latDownDistance = latDownTemp;
                latDownPoint = *it1;
            }

            centerLatIn = centerLatIn + *it1;
            contLat++;
        }
        else
        {
            medDownTemp = downMedial.getDistanceFromPoint(*it1);
            medialTemp = medialLine.getDistanceFromPoint(*it1);

            if (medDownTemp > medTopDistance)
            {
                medTopDistance = medDownTemp;
                medTopPoint = *it1;
            }

            if (medialTemp > medialDistance)
            {
                medialDistance = medialTemp;
                medPoint = *it1;
            }

            if (medDownTemp < medDownDistance || medDownDistance < 0)
            {
                medDownDistance = medDownTemp;
                medDownPoint = *it1;
            }

            centerMedIn = centerMedIn + *it1;
            contMed++;
        }
    }

    if (contLat == 0 || contMed == 0)
    {
        throw ImplantExceptionCode::CAN_NOT_DEFINE_DIFFERENT_BORDER_FOR_LATERAL_AND_MEDIAL_SIDE;
    }

    centerLatIn = centerLatIn / contLat;
    centerMedIn = centerMedIn / contMed;

    Plane planeDown;
    planeDown.init(vectorAP, downPoint);
    planeDown.reverseByPoint(centerLatIn, false);

    Point vectorLat = myMidPlane.getNormalVector();
    Point vectorDown = planeDown.getNormalVector();

    int outsideCornerPosLat = ImplantTools::GetCornerPointOnContour(points, centerLatIn, (centerLatIn + 1000000 * vectorLat), (centerLatIn + 1000000 * vectorDown), Plane(), 1, 2);
    int outsideCornerPosMed = ImplantTools::GetCornerPointOnContour(points, centerMedIn, (centerMedIn - 1000000 * vectorLat), (centerMedIn + 1000000 * vectorDown), Plane(), 1, 2);

    planeDown.movePlane(latDownPoint);
    Point latExt = planeDown.getProjectionPoint(points[outsideCornerPosLat]);

    planeDown.movePlane(medDownPoint);
    Point medExt = planeDown.getProjectionPoint(points[outsideCornerPosMed]);

    std::vector<Point> lateralOut, medialOut;
    Point centerPoint;
    separateMedialAndLateralPoints(points, midPlane, latPoint, medPoint, lateralOut, medialOut, centerPoint);

    if (lateralOut.size() == 0 || medialOut.size() == 0)
    {
        throw ImplantExceptionCode::CAN_NOT_DEFINE_DIFFERENT_BORDER_FOR_LATERAL_AND_MEDIAL_SIDE;
    }

    Plane sagitalPlane;
    sagitalPlane.init(midPlane.getNormalVector(), centerPoint);
    std::vector<Point> concaveLat, concaveMed, mainLatPoints, mainMedPoints;

    lateralOut.push_back(latExt);
    medialOut.push_back(medExt);

    ConvexHull latHull(lateralOut, pRotation);
    ConvexHull medHull(medialOut, pRotation);
    concaveLat = latHull.GetConvexHull();
    concaveMed = medHull.GetConvexHull();

    //auto tContour = ImplantTools::getContours(knee.GetFemurPoly(), currentPlane.getNormalVector(), currentPlane.getPoint());
    //ImplantTools::show(tContour, concaveLat, false);
    //ImplantTools::show(tContour, concaveMed, false);

    if (concaveLat.size() == 0)
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_CONVEX_HULL_ON_LATERAL_SIDE;
    }

    if (concaveMed.size() == 0)
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_CONVEX_HULL_ON_MEDIAL_SIDE;
    }

    bool latChange, medChange;

    std::pair<int, int> posLat, posMed;

    latChange = getConcaveMainPoints(concaveLat, currentPlane, sagitalPlane, downPoint, distanceTop, 1, mainLatPoints, posLat, 0.25);
    medChange = getConcaveMainPoints(concaveMed, currentPlane, sagitalPlane, downPoint, distanceTop, 1, mainMedPoints, posMed, 0.25);

    //tContour = ImplantTools::getContours(knee.GetFemurPoly(), currentPlane.getNormalVector(), currentPlane.getPoint());
    //ImplantTools::show(tContour, mainLatPoints, false);
    //ImplantTools::show(tContour, mainMedPoints, false);

    bool isLatBegin = true;

    if (mainLatPoints.size() > 1)
    {
        //mainLatPoints.insert(mainLatPoints.begin(), extremeDownPlane.getProjectionPoint(mainLatPoints[0]));

        if (latChange == true)
        {
            std::reverse(mainLatPoints.begin(), mainLatPoints.end());
            isLatBegin = false;

            posLat.first = mainLatPoints.size() - posLat.first - 1;
            posLat.second = mainLatPoints.size() - posLat.second - 1;
        }
    }
    else
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_MAIN_POINTS_IN_CONVEX_HULL_ON_LATERAL_SIDE;
    }

    if (mainMedPoints.size() > 1)
    {
        //mainMedPoints.insert(mainMedPoints.begin(), extremeDownPlane.getProjectionPoint(mainMedPoints[0]));

        if (medChange == true)
        {
            std::reverse(mainMedPoints.begin(), mainMedPoints.end());

            posMed.first = mainMedPoints.size() - posMed.first - 1;
            posMed.second = mainMedPoints.size() - posMed.second - 1;
        }
    }
    else
    {
        throw ImplantExceptionCode::CAN_NOT_DETERMINE_MAIN_POINTS_IN_CONVEX_HULL_ON_MEDIAL_SIDE;
    }

    std::vector<Point> allVertices;

    if (isLatBegin == true)
    {
        for (int i = 0; i < mainLatPoints.size(); i++)
        {
            allVertices.push_back(mainLatPoints[i]);
        }

        posMed.first = allVertices.size() + posMed.first;
        posMed.second = allVertices.size() + posMed.second;

        for (int i = 0; i < mainMedPoints.size(); i++)
        {
            allVertices.push_back(mainMedPoints[i]);
        }
    }
    else
    {
        for (int i = 0; i < mainMedPoints.size(); i++)
        {
            allVertices.push_back(mainMedPoints[i]);
        }

        posLat.first = allVertices.size() + posLat.first;
        posLat.second = allVertices.size() + posLat.second;

        for (int i = 0; i < mainLatPoints.size(); i++)
        {
            allVertices.push_back(mainLatPoints[i]);
        }
    }

    /*auto tContour = ImplantTools::getContours(knee.GetFemurPoly(), currentPlane.getNormalVector(), currentPlane.getPoint());
    ImplantTools::show(tContour, allVertices, true);*/

    int posInitOut, posEndOut, posInitIn, posEndIn;

    if (posLat.first < posMed.first)
    {
        posInitIn = posLat.first;
        posEndIn = posMed.first;
    }
    else
    {
        posInitIn = posMed.first;
        posEndIn = posLat.first;
    }

    if (posLat.second < posMed.second)
    {
        posInitOut = posLat.second;
        posEndOut = posMed.second;
    }
    else
    {
        posInitOut = posMed.second;
        posEndOut = posLat.second;
    }

    Point computePoint;

    for (int i = 0; i < allVertices.size() - 1; i++)
    {
        if (i >= posInitOut && i <= posEndOut)
        {
            vertices.push_back(allVertices[i]);
            continue;
        }

        if (distanceSide != 0)
        {
            computePoint = movePointAtNormal(allVertices[i], allVertices[i + 1], pRotation, distanceSide);
            vertices.push_back(computePoint);

            if (i == allVertices.size() - 2)
            {
                computePoint = movePointAtNormal(allVertices[i], allVertices[i + 1], pRotation, distanceSide, true);
                vertices.push_back(computePoint);
            }
        }
        else
        {
            vertices.push_back(allVertices[i]);
        }
    }

    vertices = ConvexHull::interpolateSpline(vertices, amount);
}

void FemurImplantMatch::getVerticesCDE(const std::vector<Point>& points, const Point& downPoint, const Point& lateralPoint, const Point& medialPoint, const Point& topPoint, const Point& centerP1, const Point& centerP2, double distanceSide, double distanceTop, double angleLat, double angleMed, std::vector<Point>& vertices) const
{
    Point directVectorAP = centerP2 - centerP1;
    Line lateralLine(directVectorAP, lateralPoint);
    Line medialLine(directVectorAP, medialPoint);

    if (lateralLine.isPointBelongToLine(medialPoint))
    {
        throw ImplantExceptionCode::CAN_NOT_DEFINE_DIFFERENT_BORDER_FOR_LATERAL_AND_MEDIAL_SIDE;
    }

    Line topLineTemp = lateralLine.getPerpendicularLine(medialPoint);
    Line topLine(topLineTemp.getDirectVector(), topPoint);

    if (topLine.isPointBelongToLine(downPoint))
    {
        throw ImplantExceptionCode::CAN_NOT_DEFINE_DIFFERENT_BORDER_FOR_INFERIOR_AND_SUPERIOR_SIDE;
    }

    directVectorAP = topLine.getProjectPoint(downPoint) - downPoint;
    directVectorAP.normalice();

    Plane transversePlane;
    transversePlane.init(directVectorAP, topPoint);
    Point lateralIntercept = transversePlane.getInterceptionLinePoint(lateralLine);
    Point medialIntercept = transversePlane.getInterceptionLinePoint(medialLine);
    Point transverseVector = lateralLine.getProjectPoint(medialPoint) - medialPoint;
    transverseVector.normalice();

    Point latObliquePoint = lateralIntercept + transverseVector - (tan(angleLat) * directVectorAP);
    Point medObliquePoint = medialIntercept - transverseVector - (tan(angleMed) * directVectorAP);

    lateralLine.setDirectVector((lateralIntercept - latObliquePoint));
    medialLine.setDirectVector((medialIntercept - medObliquePoint));

    double lateralDistance = -1;
    double medialDistance = -1;
    double topDistance = -1;
    double downDistance = -1;
    double lateralTemp, medialTemp, toTemp, downTemp;
    Point sideLateralPoint, sideMedialPoint, sideTopPoint, midPoint, sideDownPoint;

    std::vector<Point>::const_iterator it1, it2;
    it1 = points.begin();
    it2 = points.end();
    for (; it1 != it2; ++it1)
    {
        lateralTemp = lateralLine.getDistanceFromPoint(*it1);
        medialTemp = medialLine.getDistanceFromPoint(*it1);
        toTemp = topLine.getDistanceFromPoint(*it1);
        downTemp = topLine.getDistanceFromPoint(*it1);

        if (lateralTemp < lateralDistance || lateralDistance < 0)
        {
            lateralDistance = lateralTemp;
            sideLateralPoint = *it1;
        }

        if (medialTemp < medialDistance || medialDistance < 0)
        {
            medialDistance = medialTemp;
            sideMedialPoint = *it1;
        }

        if (toTemp < topDistance || topDistance < 0)
        {
            topDistance = toTemp;
            sideTopPoint = *it1;
        }

        if (downTemp > downDistance)
        {
            downDistance = downTemp;
            sideDownPoint = *it1;
        }

        midPoint = midPoint + *it1;
    }

    Point latPerpendicularV = lateralLine.getProjectPoint(sideLateralPoint) - sideLateralPoint;
    Point medPerpendicularV = medialLine.getProjectPoint(sideMedialPoint) - sideMedialPoint;
    Point topPerpendicularV = topLine.getProjectPoint(sideTopPoint) - sideTopPoint;
    latPerpendicularV.normalice();
    medPerpendicularV.normalice();
    topPerpendicularV.normalice();

    lateralLine.setPoint(sideLateralPoint + distanceSide * latPerpendicularV);
    medialLine.setPoint(sideMedialPoint + distanceSide * medPerpendicularV);
    Point finalTopPoint = sideTopPoint + distanceTop * topPerpendicularV;

    std::vector<Point> verticesTemp;

    transversePlane.movePlane(sideDownPoint);

    Point downLat = transversePlane.getInterceptionLinePoint(lateralLine);
    Point downMed = transversePlane.getInterceptionLinePoint(medialLine);

    transversePlane.movePlane(finalTopPoint);

    Point topLat = transversePlane.getInterceptionLinePoint(lateralLine);
    Point topMed = transversePlane.getInterceptionLinePoint(medialLine);

    vertices = { downLat, topLat, topMed, downMed };
}

/*
void FemurImplantMatch::getVerticesCDE(const std::vector<Point>& points, const Point& downPoint, const Point& lateralPoint, const Point& medialPoint, const Point& topPoint, const Point& centerP1, const Point& centerP2, std::vector<Point>& vertices, double distance) const
{
    double tan75 = 3.7;
    Point directVectorAP = centerP2 - centerP1;
    Line lateralLine(directVectorAP, lateralPoint);
    Line medialLine(directVectorAP, medialPoint);

    if (lateralLine.isPointBelongToLine(medialPoint))
    {
        throw ImplantExceptionCode::CAN_NOT_DEFINE_DIFFERENT_BORDER_FOR_LATERAL_AND_MEDIAL_SIDE;
    }

    Line topLineTemp = lateralLine.getPerpendicularLine(medialPoint);
    Line topLine(topLineTemp.getDirectVector(), topPoint);

    if (topLine.isPointBelongToLine(downPoint))
    {
        throw ImplantExceptionCode::CAN_NOT_DEFINE_DIFFERENT_BORDER_FOR_INFERIOR_AND_SUPERIOR_SIDE;
    }

    directVectorAP = topLine.getProjectPoint(downPoint) - downPoint;
    directVectorAP.normalice();

    Plane transversePlane;
    transversePlane.init(directVectorAP, topPoint);
    Point lateralIntercept = transversePlane.getInterceptionLinePoint(lateralLine);
    Point medialIntercept = transversePlane.getInterceptionLinePoint(medialLine);
    Point transverseVector = lateralLine.getProjectPoint(medialPoint) - medialPoint;
    transverseVector.normalice();

    Point topLatObliquePoint = lateralIntercept - transverseVector;
    Point topMedObliquePoint = medialIntercept + transverseVector;

    Point latObliquePoint = lateralIntercept - directVectorAP * tan75;
    Point medObliquePoint = medialIntercept - directVectorAP * tan75;

    lateralLine.setDirectVector((topLatObliquePoint - latObliquePoint));
    medialLine.setDirectVector((topMedObliquePoint - medObliquePoint));

    double myEpsilo = 0.01;
    double lateralDistance = -1.0;
    double medialDistance = -1.0;
    double topDistance = -1.0;
    double downDistance = -1.0;
    double lateralTemp, medialTemp, toTemp, downTemp;
    Point sideLateralPoint, sideMedialPoint, sideTopPoint, midPoint, sideDownPoint;

    std::vector<Point>::const_iterator it1, it2;
    it1 = points.begin();
    it2 = points.end();
    for (; it1 != it2; ++it1)
    {
        lateralTemp = 1.0 / (lateralLine.getDistanceFromPoint(*it1) + myEpsilo);
        medialTemp = 1.0 / (medialLine.getDistanceFromPoint(*it1) + myEpsilo);
        toTemp = 1.0 / (topLine.getDistanceFromPoint(*it1) + myEpsilo);
        downTemp = topLine.getDistanceFromPoint(*it1);
        if (lateralTemp > lateralDistance)
        {
            lateralDistance = lateralTemp;
            sideLateralPoint = *it1;
        }

        if (medialTemp > medialDistance)
        {
            medialDistance = medialTemp;
            sideMedialPoint = *it1;
        }

        if (toTemp > topDistance)
        {
            topDistance = toTemp;
            sideTopPoint = *it1;
        }

        if (downTemp > downDistance)
        {
            downDistance = downTemp;
            sideDownPoint = *it1;
        }

        midPoint = midPoint + *it1;
    }

    Point latPerpendicularV = lateralLine.getProjectPoint(sideLateralPoint) - sideLateralPoint;
    Point medPerpendicularV = medialLine.getProjectPoint(sideMedialPoint) - sideMedialPoint;
    Point topPerpendicularV = topLine.getProjectPoint(sideTopPoint) - sideTopPoint;
    latPerpendicularV.normalice();
    medPerpendicularV.normalice();
    topPerpendicularV.normalice();

    lateralLine.setPoint(sideLateralPoint + distance * latPerpendicularV);
    medialLine.setPoint(sideMedialPoint + distance * medPerpendicularV);
    Point finalTopPoint = sideTopPoint + distance * topPerpendicularV;

    std::vector<Point> verticesTemp;

    transversePlane.movePlane(sideDownPoint);

    Point downLat = transversePlane.getInterceptionLinePoint(lateralLine);
    Point downMed = transversePlane.getInterceptionLinePoint(medialLine);

    transversePlane.movePlane(finalTopPoint);

    Point topLat = transversePlane.getInterceptionLinePoint(lateralLine);
    Point topMed = transversePlane.getInterceptionLinePoint(medialLine);

    vertices = { downLat, topLat, topMed, downMed };
}
*/

bool FemurImplantMatch::getConcaveMainPoints(const std::vector<Point>& points, const Plane& currentPlane, const Plane& sagitalPlane, const Point& downPoint, int pDistanceTop, int pDistanceSideIn, std::vector<Point>& result, std::pair<int, int>& pInOutPos, double pMoveOnCenter) const
{
    Point extremeOutSide, extremeInSide, centerPoint, extremeOutSideTemp, extremeTop;
    double distance, distanceTempInside, distanceTempOutside, distanceTempTop;
    distanceTempInside = 9999999999.0;
    distanceTempOutside = -1;
    distanceTempTop = -1;
    bool changeDirection = false;
    auto it1 = points.begin();
    auto it2 = points.end();

    Point vectorAP = currentPlane.getNormalVector().cross(sagitalPlane.getNormalVector());
    Plane downPlane;
    downPlane.init(vectorAP, downPoint);

    for (; it1 != it2; ++it1)
    {
        centerPoint = centerPoint + *it1;
        distance = sagitalPlane.getDistanceFromPoint(*it1);
        if (distance < distanceTempInside)
        {
            distanceTempInside = distance;
            extremeInSide = *it1;
        }

        if (distance > distanceTempOutside)
        {
            distanceTempOutside = distance;
            extremeOutSide = *it1;
        }

        distance = downPlane.getDistanceFromPoint(*it1);

        if (distance > distanceTempTop)
        {
            distanceTempTop = distance;
            extremeTop = *it1;
        }
    }

    downPlane.movePlane(extremeTop);
    Line insideLine(downPlane.getNormalVector(), extremeInSide);
    Point intercepTop = downPlane.getInterceptionLinePoint(insideLine);
    downPlane.movePlane(downPoint);

    centerPoint = centerPoint / double(points.size());

    vectorAP = downPlane.getProjectionPoint(centerPoint) - centerPoint;
    vectorAP.normalice();

    Point vectorTrans = extremeOutSide - sagitalPlane.getProjectionPoint(extremeOutSide);
    vectorTrans.normalice();

    Point extremeA = centerPoint + 100 * vectorTrans;
    Point extremeB = centerPoint + 100 * vectorAP;

    int outsideCornerPos = ImplantTools::GetCornerPointOnContour(points, centerPoint, extremeA, extremeB);

    if (outsideCornerPos < 0)
    {
        throw ImplantExceptionCode::CAN_NOT_DETECT_REFERENCE_POINT_AT_INFERIOR_CONTOUR_CORNER;
    }

    std::vector<Point> tempPoints = points;

    std::rotate(tempPoints.begin(), tempPoints.begin() + outsideCornerPos, tempPoints.end());
    Plane oblique = currentPlane.getPerpendicularPlane(centerPoint, points[outsideCornerPos]);
    oblique.reverseByPoint(extremeA);
    if (oblique.eval(tempPoints[2]) < 0)
    {
        std::reverse(tempPoints.begin(), tempPoints.end());
        changeDirection = true;
    }

    Plane transversal = currentPlane.getPerpendicularPlane(centerPoint, extremeTop);
    transversal.reverseByPoint(extremeInSide, false);

    for (int i = 1; i < tempPoints.size() - 1; i++) // (tempPoints.size() - 1) because the last point can be so close to the plane that its evaluation can be positive or negative. It is not a safe point.
    {
        if (i == 1)
        {
            result.push_back(tempPoints[0]);
        }

        if (oblique.eval(tempPoints[i]) > 0 && transversal.eval(tempPoints[i]) > 0)
        {
            result.push_back(tempPoints[i]);
        }
    }

    //auto tContour = ImplantTools::getContours(knee.GetFemurPoly(), currentPlane.getNormalVector(), currentPlane.getPoint());
    //ImplantTools::show(tContour, result, false);

    Point lastTopPoint = result[result.size() - 1];

    double slopeAngle = 92;
    downPlane.movePlane(extremeTop);
    downPlane.reverseByPoint(downPoint, false);
    downPlane.movePlaneOnNormal(pDistanceTop);
    int pDataSidePos;
    Line myLineInside = ImplantTools::GetSquareCornerFeatures(lastTopPoint, intercepTop, insideLine.getProjectPoint(centerPoint), points, pDataSidePos, slopeAngle);
    Point intercepTopIn = downPlane.getInterceptionLinePoint(myLineInside);
    Point sideIn = myLineInside.getProjectPoint(centerPoint);

    sideIn = sideIn + pMoveOnCenter * (intercepTopIn - sideIn);

    insideLine.setPoint(extremeOutSide);
    Line myLineOutSide = ImplantTools::GetSquareCornerFeatures(lastTopPoint, downPlane.getInterceptionLinePoint(insideLine), extremeOutSide, result, pDataSidePos, slopeAngle);
    Point intercepTopOut = downPlane.getInterceptionLinePoint(myLineOutSide);
    Point sideOut = result[pDataSidePos];

    while (result.size() > pDataSidePos + 1)
    {
        result.pop_back();
    }

    //tContour = ImplantTools::getContours(knee.GetFemurPoly(), currentPlane.getNormalVector(), currentPlane.getPoint());
    //ImplantTools::show(tContour, result, false);

    for (double i = 0.2; i < 0.9; i += 0.2)
    {
        Point temp = sideOut + i * (intercepTopOut - sideOut);
        result.push_back(temp);
    }
    pInOutPos.second = result.size();
    result.push_back(intercepTopOut);

    for (double i = 0.2; i < 0.9; i += 0.2)
    {
        Point temp = intercepTopOut + i * (intercepTopIn - intercepTopOut);
        result.push_back(temp);
    }

    pInOutPos.first = result.size();
    result.push_back(intercepTopIn);

    Point vectorInShift = intercepTopIn - intercepTopOut;
    vectorInShift.normalice();

    for (double i = 0.2; i < 0.9; i += 0.2)
    {
        Point temp = intercepTopIn + i * (sideIn - intercepTopIn) + (pDistanceSideIn * vectorInShift);
        result.push_back(temp);
    }

    result.push_back((sideIn + pDistanceSideIn * vectorInShift));

    //tContour = ImplantTools::getContours(knee.GetFemurPoly(), currentPlane.getNormalVector(), currentPlane.getPoint());
    //ImplantTools::show(tContour, result, false);

    //tContour = ImplantTools::getContours(knee.GetFemurPoly(), currentPlane.getNormalVector(), currentPlane.getPoint());
    //ImplantTools::show(tContour, result, true);

    /*
    if (intercepTop != lastTopPoint)
    {
        for (double i = 0.2; i < 0.9; i += 0.2)
        {
            Point temp = lastTopPoint + i * (intercepTop - lastTopPoint);
            result.push_back(temp);
        }
        result.push_back(intercepTop);
    }

    tContour = ImplantTools::getContours(knee.GetFemurPoly(), currentPlane.getNormalVector(), currentPlane.getPoint());
    ImplantTools::show(tContour, result, false);

    Point centerProj = insideLine.getProjectPoint(centerPoint);

    for (double i = 0.2; i < 0.9; i += 0.2)
    {
        Point temp = intercepTop + i * (centerProj - intercepTop);
        result.push_back(temp);
    }
    result.push_back(centerProj);

    tContour = ImplantTools::getContours(knee.GetFemurPoly(), currentPlane.getNormalVector(), currentPlane.getPoint());
    ImplantTools::show(tContour, result, false);

    downPlane.reverseByPoint(downPoint, false);
    Point supPoint = centerPoint + 1000 * downPlane.getNormalVector();

    ImplantTools::squareCorner(supPoint, centerPoint, extremeOutSide, result);

    tContour = ImplantTools::getContours(knee.GetFemurPoly(), currentPlane.getNormalVector(), currentPlane.getPoint());
    ImplantTools::show(tContour, result, false);
    */

    return changeDirection;
}

std::pair<int, int> FemurImplantMatch::getSupExternalCornerPos(const std::vector<Point>& allPoints, const std::vector<Point>& contour, const Plane& currentPlane, const Plane& sagitalPlane, const Point& downPoint) const
{
    Point extremeOutSide, extremeInSide, centerPoint;
    double distance, distanceTempInside, distanceTempOutside;
    distanceTempInside = 9999999.0;
    distanceTempOutside = -1;

    bool changeDirection = false;
    auto it1 = contour.begin();
    auto it2 = contour.end();

    for (; it1 != it2; ++it1)
    {
        centerPoint = centerPoint + *it1;
        distance = sagitalPlane.getDistanceFromPoint(*it1);
        if (distance < distanceTempInside)
        {
            distanceTempInside = distance;
            extremeInSide = *it1;
        }

        if (distance > distanceTempOutside)
        {
            distanceTempOutside = distance;
            extremeOutSide = *it1;
        }
    }
    centerPoint = centerPoint / double(contour.size());

    Point vectorAP = currentPlane.getNormalVector().cross(sagitalPlane.getNormalVector());
    Plane downPlane;
    downPlane.init(vectorAP, downPoint);
    vectorAP = centerPoint - downPlane.getProjectionPoint(centerPoint);
    vectorAP.normalice();

    Point vectorTrans = extremeOutSide - sagitalPlane.getProjectionPoint(extremeOutSide);
    vectorTrans.normalice();

    Point extremeOut = centerPoint + 100 * vectorTrans;
    Point extremeTop = centerPoint + 100 * vectorAP;
    Point extremeIn = centerPoint - 100 * vectorTrans;

    Plane helpPlane = sagitalPlane;
    helpPlane.reverseByPoint(extremeOut);

    int outsideCornerPos = ImplantTools::GetCornerPointOnContour(allPoints, centerPoint, extremeOut, extremeTop);
    int insideCornerPos = ImplantTools::GetCornerPointOnContour(allPoints, centerPoint, extremeIn, extremeTop, helpPlane);

    if (outsideCornerPos < 0 || insideCornerPos < 0)
    {
        throw ImplantExceptionCode::CAN_NOT_DETECT_REFERENCE_POINT_AT_SUPERIOR_CONTOUR_CORNER;
    }
    return std::make_pair(insideCornerPos, outsideCornerPos);
}

void FemurImplantMatch::separateMedialAndLateralPoints(const std::vector<Point>& points, const Plane& sagitalPlane, const Point& extremeLatPoint, const Point& extremeMedPoint, std::vector<Point>& lateralOut, std::vector<Point>& medialOut, Point& centerPointOut) const
{
    Plane myMovePlane;
    myMovePlane.init(sagitalPlane.getNormalVector(), extremeLatPoint);
    myMovePlane.normalizeNormalVector();

    Point moveVector = extremeMedPoint - myMovePlane.getProjectionPoint(extremeMedPoint);
    moveVector.normalice();

    if (myMovePlane.eval(extremeMedPoint) > 0)
    {
        myMovePlane.reverse();
    }

    std::vector<Point> latTempPoints, medTempPoints;
    Point movePoint = extremeLatPoint;

    double latDistance, medDistance, distanceTemp;

    Point nearLat, nearMed, temp;

    PointsBorder edge;

    double maxDistance = 9999999.0;

    do
    {
        latDistance = maxDistance;
        medDistance = maxDistance;
        latTempPoints.clear();
        medTempPoints.clear();
        movePoint = movePoint + 5.0 * moveVector;
        myMovePlane.movePlane(movePoint);
        auto it1 = points.begin();
        auto it2 = points.end();
        for (; it1 != it2; ++it1)
        {
            temp = *it1;
            distanceTemp = myMovePlane.eval(temp);
            if (distanceTemp >= 0)
            {
                latTempPoints.push_back(temp);
                if (abs(distanceTemp) < latDistance)
                {
                    latDistance = abs(distanceTemp);
                    nearLat = temp;
                }
            }
            else if (distanceTemp < 0)
            {
                medTempPoints.push_back(temp);
                if (abs(distanceTemp) < medDistance)
                {
                    medDistance = abs(distanceTemp);
                    nearMed = temp;
                }
            }
        }

        if (latDistance < maxDistance && medDistance < maxDistance)
        {
            double edgeDistance = edge.latDistance + edge.medDistance;
            if (edgeDistance < latDistance + medDistance)
            {
                edge.latDistance = latDistance;
                edge.medDistance = medDistance;
                edge.latNearPoint = nearLat;
                edge.medNearPoint = nearMed;
                edge.latPoints = latTempPoints;
                edge.medPoints = medTempPoints;
            }
        }

    } while (latDistance < maxDistance && medDistance < maxDistance);

    if (edge.latDistance > 0 && edge.medDistance > 0 && (edge.latDistance + edge.medDistance) > 1.0)
    {
        myMovePlane.movePlane(edge.latNearPoint);
        centerPointOut = (myMovePlane.getProjectionPoint(edge.medNearPoint) + edge.medNearPoint) / 2.0;
        lateralOut = edge.latPoints;
        medialOut = edge.medPoints;
    }
    else
    {
        //std::cout << "Distance lat: " << edge.latDistance << "Distance med: " << edge.medDistance <<std::endl;
        throw ImplantExceptionCode::CAN_NOT_DEFINE_DIFFERENT_BORDER_FOR_LATERAL_AND_MEDIAL_SIDE;
    }

}

std::vector<PointTypeITK> FemurImplantMatch::increaseVectorToAmount(const std::vector<Point>& points, int amount) const
{
    std::vector<PointTypeITK> result;
    if (points.size() >= amount || points.size() <= 1 || amount < 3)
    {
        auto it1 = points.begin();
        auto it2 = points.end();

        for (; it1 != it2; ++it1)
        {
            result.push_back((*it1).ToITKPoint());
        }
        return result;
    }

    /*int interPoints = amount / (points.size() - 1);
    int lastPos = 0;
    int stillPoints = amount - points.size();
    int interAmount = stillPoints / interPoints;
    int rest = stillPoints % interPoints;
    Point a, b, c;
    double coef = 1.0 / double(interPoints + 1);
    for (int i = 0; i < points.size() - 1; i++)
    {
        lastPos = i;
        a = points[i];
        b = points[i + 1];
        result.push_back(a.ToITKPoint());
        for (int j = 1; j <= interPoints; j++)
        {
            c = a + double(j) * coef * (b - a);
            result.push_back(c.ToITKPoint());
        }
        interAmount--;
        if (interAmount == 0)
        {
            break;
        }
    }*/

    int intervals = points.size() - 2;
    if (intervals == 0)
    {
        intervals = 1;
    }
    int stillPoints = amount - points.size();
    int interAmount = stillPoints / intervals;
    int lastPos = 0;
    int rest = stillPoints % intervals;
    Point a, b, c;
    double coef = 1.0 / double(interAmount + 1);
    double coefRest = 1.0 / double(rest + 1);
    for (int i = 0; i < points.size() - 1; i++)
    {
        lastPos = i;
        a = points[i];
        b = points[i + 1];
        result.push_back(a.ToITKPoint());
        for (int j = 1; j <= interAmount; j++)
        {
            c = a + double(j) * coef * (b - a);
            result.push_back(c.ToITKPoint());
        }
        intervals--;
        if (intervals == 0)
        {
            break;
        }
    }

    for (int i = lastPos + 1; i < points.size() - 1; i++)
    {
        if (i == lastPos + 1)
        {
            a = points[i];
            b = points[i + 1];

            result.push_back(a.ToITKPoint());
            for (int j = 1; j <= rest; j++)
            {
                c = a + double(j) * coefRest * (b - a);
                result.push_back(c.ToITKPoint());
            }
        }
        else
        {
            result.push_back((points[i]).ToITKPoint());
        }
    }
    result.push_back((points[points.size() - 1]).ToITKPoint());
    return result;
}

Point FemurImplantMatch::movePointAtNormal(const Point& movePoint, const Point& nextPoint, const cv::Mat& rotationZ, double distance, bool changeMove, bool clockWise) const
{
    cv::Mat moveMat = rotationZ * movePoint.ToMatPoint();
    cv::Mat nextMat = rotationZ * nextPoint.ToMatPoint();

    Point moveP = Point(moveMat);
    Point nextP = Point(nextMat);
    Point vector;
    if (clockWise == true)
    {
        vector = nextP - moveP;
    }
    else
    {
        vector = moveP - nextP;
    }

    cv::Point2d perpendicular2d(vector.y, -vector.x);
    perpendicular2d = perpendicular2d / sqrt(perpendicular2d.dot(perpendicular2d));
    Point perpendicular(perpendicular2d.x, perpendicular2d.y, 0);
    Point finalMove;
    if (changeMove == true)
    {
        finalMove = nextP + distance * perpendicular;
    }
    else
    {
        finalMove = moveP + distance * perpendicular;
    }

    cv::Mat resultMat = rotationZ.inv() * finalMove.ToMatPoint();
    return Point(resultMat);
}