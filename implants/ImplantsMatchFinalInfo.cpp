#include "ImplantsMatchFinalInfo.hpp"
#include "LegAngle.hpp"
#include "BoneRotulaPath.hpp"
#include "ImplantTools.hpp"
#include <itkVersorRigid3DTransform.h>

using namespace TKA::IMPLANTS;

ImplantsMatchFinalInfo::ImplantsMatchFinalInfo(Knee* pKnee, const FemurImplant pFemurImplant, const TibiaImplant pTibiaImplant, const itk::Rigid3DTransform<>::Pointer pImplantToBoneFemurTransform, const itk::Rigid3DTransform<>::Pointer pImplantToBoneTibiaTransform)
{
    knee = pKnee;
    femurRotation = Rigid3DTransformToCVRotation(pImplantToBoneFemurTransform);
    femurTranslation = Rigid3DTransformToCVTranslation(pImplantToBoneFemurTransform);
    tibiaRotation = Rigid3DTransformToCVRotation(pImplantToBoneTibiaTransform);
    tibiaTranslation = Rigid3DTransformToCVTranslation(pImplantToBoneTibiaTransform);
    femurImplant = pFemurImplant;
    tibiaImplant = pTibiaImplant;
    femurCartilage = pKnee->getFemurCartilage();
    tibiaCartilage = pKnee->getTibiaCartilage();
    /*femurLatThickness = knee->getImplantInfo().femurVarusLateralThickness;
    femurMedialThickness = knee->getImplantInfo().femurValgusMedialThickness;
    tibiaLatThickness = knee->getImplantInfo().tibiaVarusLateralThickness;
    tibiaMedialThickness = knee->getImplantInfo().tibiaValgusMedialThickness;*/
    //femurCoordenate = new CoordenateSystemFemur(pKnee->getHipCenter(), pKnee->getLateralEpicondyle(), pKnee->getMedialEpicondylePerp(), pKnee->getFemurKneeCenter());

    //BoneRotulaPath object(pKnee);
    this->boneKneeCapPath = pKnee->getKneeGroovePath();
    this->implantKneeCapPath = pFemurImplant.GetKneeCapPath();

    updateFemurImplantVectors();
    updateTibiaImplantVectors();
}

void ImplantsMatchFinalInfo::updateFemurImplantVectors()
{
    cv::Mat axisTemp = femurRotation * (femurImplant.getDirectVectorFemurAxis().ToMatPoint());
    cv::Mat apTemp = femurRotation * (femurImplant.getDirectVectorAP().ToMatPoint());
    cv::Mat TEATemp = femurRotation * (femurImplant.getDirectVectorTEA().ToMatPoint());
    cv::Mat middleTemp = femurRotation * (femurImplant.getMidPlaneInterceptionPoint().ToMatPoint()) + femurTranslation;

    femurAxisImplantVector = Point(axisTemp);
    femurAPImplantVector = Point(apTemp);
    femurTEAImplantVector = Point(TEATemp);

    femurAxisImplantVector.normalice();
    femurAPImplantVector.normalice();
    femurTEAImplantVector.normalice();
    femurImplantMiddlePoint = Point(middleTemp);
}

void ImplantsMatchFinalInfo::updateTibiaImplantVectors()
{
    cv::Mat axisTemp = tibiaRotation * (tibiaImplant.getTibiaNormalVector().ToMatPoint());
    cv::Mat apTemp = tibiaRotation * (tibiaImplant.getTibiaVectorAP().ToMatPoint());
    cv::Mat TEATemp = tibiaRotation * (tibiaImplant.getTibiaVectorTEA().ToMatPoint());
    cv::Mat middleTemp = tibiaRotation * (tibiaImplant.getCentralPoint().ToMatPoint()) + tibiaTranslation;

    tibiaAxisImplantVector = Point(axisTemp);
    tibiaAPImplantVector = Point(apTemp);
    tibiaTEAImplantVector = Point(TEATemp);

    tibiaAxisImplantVector.normalice();
    tibiaAPImplantVector.normalice();
    tibiaTEAImplantVector.normalice();
    tibiaImplantMiddlePoint = Point(middleTemp);
}

ImplantsMatchFinalInfo::~ImplantsMatchFinalInfo()
{
    /*delete femurCoordenate;
    femurCoordenate = NULL;*/
}

const itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::setFemurImplant(const FemurImplant& pFemurImplant)
{
    femurImplant = pFemurImplant;
    this->implantKneeCapPath = pFemurImplant.GetKneeCapPath();

    std::vector<cv::Point3d> implantVectors;
    std::vector<cv::Point3d> kneeVectors;

    Point implantTEA = femurImplant.getDirectVectorTEA();
    Point implantFemurAxis = femurImplant.getDirectVectorFemurAxis();
    Point implantAP = femurImplant.getDirectVectorAP();
    implantVectors.push_back(implantTEA.ToCVPoint());
    implantVectors.push_back(implantFemurAxis.ToCVPoint());
    implantVectors.push_back(implantAP.ToCVPoint());

    kneeVectors.push_back(femurTEAImplantVector.ToCVPoint());
    kneeVectors.push_back(femurAxisImplantVector.ToCVPoint());
    kneeVectors.push_back(femurAPImplantVector.ToCVPoint());

    cv::Mat implantMatrix = cv::Mat(implantVectors.size(), 3, CV_64F, implantVectors.data());
    cv::Mat kneeMatrix = cv::Mat(kneeVectors.size(), 3, CV_64F, kneeVectors.data());

    cv::Mat inverse = (implantMatrix.t()).inv();

    femurRotation = (kneeMatrix.t()) * inverse;
    femurTranslation = femurImplantMiddlePoint.ToMatPoint() - (femurRotation * femurImplant.getMidPlaneInterceptionPoint().ToMatPoint());

    return getITKFemurTransform();
}

const itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::setTibiaImplant(const TibiaImplant& pTibiaImplant)
{
    tibiaImplant = pTibiaImplant;

    std::vector<cv::Point3d> implantVectors;
    std::vector<cv::Point3d> kneeVectors;

    Point implantNormal = tibiaImplant.getTibiaNormalVector();
    Point implantAP = tibiaImplant.getTibiaVectorAP();
    Point implantCross = tibiaImplant.getTibiaVectorTEA();
    implantVectors.push_back(implantNormal.ToCVPoint());
    implantVectors.push_back(implantAP.ToCVPoint());
    implantVectors.push_back(implantCross.ToCVPoint());

    kneeVectors.push_back(tibiaAxisImplantVector.ToCVPoint());
    kneeVectors.push_back(tibiaAPImplantVector.ToCVPoint());
    kneeVectors.push_back(tibiaTEAImplantVector.ToCVPoint());

    cv::Mat implantMatrix = cv::Mat(implantVectors.size(), 3, CV_64F, implantVectors.data());
    cv::Mat kneeMatrix = cv::Mat(kneeVectors.size(), 3, CV_64F, kneeVectors.data());

    cv::Mat inverse = (implantMatrix.t()).inv();

    tibiaRotation = (kneeMatrix.t()) * inverse;
    tibiaTranslation = tibiaImplantMiddlePoint.ToMatPoint() - (tibiaRotation * tibiaImplant.getCentralPoint().ToMatPoint());

    return getITKTibiaTransform();
}

/*
void ImplantsMatchFinalInfo::setGeneralImplantInfo(const ImplantInfo pImplant)
{
    knee->setImplantInfo(pImplant);
    
    femurLatThickness = knee->getImplantInfo().femurVarusLateralThickness;
    femurMedialThickness = knee->getImplantInfo().femurValgusMedialThickness;
    tibiaLatThickness = knee->getImplantInfo().tibiaVarusLateralThickness;
    tibiaMedialThickness = knee->getImplantInfo().tibiaValgusMedialThickness;
    
}
*/

std::vector<PointTypeITK> ImplantsMatchFinalInfo::getBoneKneeCapPath() const
{
    std::vector<PointTypeITK> pointsOut;
    int tSize = boneKneeCapPath.size();
    for (int i = 0; i < tSize; i++)
    {
        pointsOut.push_back(boneKneeCapPath[i].ToITKPoint());
    }

    return pointsOut;
}

std::vector<PointTypeITK> ImplantsMatchFinalInfo::getImplantKneeCapPath() const
{
    std::vector<PointTypeITK> pointsOut;
    int tSize = implantKneeCapPath.size();
    for (int i = 0; i < tSize; i++)
    {
        /*cv::Mat newPoint = femurRotation * (implantKneeCapPath[i].ToMatPoint()) + femurTranslation;
        Point proj = femurCoordenate->getPointCoordenate(Point(newPoint));
        PointTypeITK itkPoint;
        itkPoint[0] = proj.x;
        itkPoint[1] = proj.y;
        itkPoint[2] = proj.z;
        pointsOut.push_back(itkPoint);*/

        cv::Mat newPoint = femurRotation * (implantKneeCapPath[i].ToMatPoint()) + femurTranslation;
        Point proj = Point(newPoint);
        pointsOut.push_back(proj.ToITKPoint());
    }

    return pointsOut;
}

const itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::setTibiaSlopeAngle(double angle)
{
    double myAngle = (angle - GetTibiaImplantSlopeAngle()) * PI / 180.0;

    Point crossVector = knee->getNormalVectorTibiaPlane().cross(knee->getTibiaDirectVectorAP());

    Plane sagital;
    sagital.init(crossVector, knee->getMedialPlateau());

    Point axisRotation = knee->getLateralPlateau() - sagital.getProjectionPoint(knee->getLateralPlateau());

    if (knee->getIsRight() == false)
    {
        axisRotation = (-1.0) * axisRotation;
    }

    cv::Mat rotation = ImplantTools::getRotateMatrix(axisRotation, myAngle);
    tibiaRotation = rotation * tibiaRotation;

    knee->setTibiaSlope(angle);

    return getITKTibiaTransform();
}

double ImplantsMatchFinalInfo::GetTibiaImplantSlopeAngle() const
{
    Point tibiaAxis = knee->getTibiaKneeCenter() - knee->getAnkleCenter();
    Plane tibiaHelp;
    tibiaHelp.init(tibiaAxis, knee->getTibiaKneeCenter());

    Point boneAP = (-1.0) * (knee->getTibiaDirectVectorAP());
    boneAP = tibiaHelp.getProjectionVector(boneAP);

    Point implantAP = (-1.0) * (tibiaImplant.getTibiaVectorAP());
    cv::Mat implantAPMat = tibiaRotation * (implantAP.ToMatPoint());
    implantAP = Point(implantAPMat);

    Plane sagital;
    sagital.init(knee->getTibiaVectorTEA(), knee->getTibiaKneeCenter());

    if (tibiaHelp.eval(knee->getAnkleCenter()) < 0)
    {
        tibiaHelp.reverse();
    }

    boneAP = sagital.getProjectionVector(boneAP);
    implantAP = sagital.getProjectionVector(implantAP);

    LegAngle legAngle;
    double angle = legAngle.getAngleBetweenVectors(boneAP, implantAP);

    if (!(angle > 0 && angle < 180))
    {
        return angle;
    }

    Point newPoint = knee->getTibiaKneeCenter() + implantAP;

    if (tibiaHelp.eval(newPoint) > 0)
    {
        return angle;
    }
    else
    {
        return -angle;
    }
}

const itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::setTibiaRotationAngle(double angle)
{
    double myAngle = (angle - GetTibiaImplantRotationAngle()) * PI / 180.0;

    Point axisRotation = knee->getNormalVectorTibiaPlane();

    if (knee->getIsRight() == true)
    {
        axisRotation = (-1.0) * axisRotation;
    }

    cv::Mat rotation = ImplantTools::getRotateMatrix(axisRotation, myAngle);
    tibiaRotation = rotation * tibiaRotation;

    return getITKTibiaTransform();
}

double ImplantsMatchFinalInfo::GetTibiaImplantRotationAngle() const
{
    Plane transversal;
    transversal.init(knee->getNormalVectorTibiaPlane(), knee->getTibiaKneeCenter());

    Point tibiaVectorTEA = knee->getNormalVectorTibiaPlane().cross(knee->getTibiaDirectVectorAP());

    Plane sagital;
    sagital.init(tibiaVectorTEA, knee->getLateralPlateau());

    if (sagital.eval(knee->getMedialPlateau()) > 0)
    {
        sagital.reverse();
    }

    Point vectorAPProj = transversal.getProjectionVector(knee->getTibiaDirectVectorAP());
    cv::Mat tibiaImplantVectorMat = tibiaRotation * (tibiaImplant.getTibiaVectorAP().ToMatPoint());
    Point tibiaImplantVector = Point(tibiaImplantVectorMat);
    tibiaImplantVector = transversal.getProjectionVector(tibiaImplantVector);
    LegAngle legAngle;

    double angle = legAngle.getAngleBetweenVectors(vectorAPProj, tibiaImplantVector);

    if (!(angle > 0 && angle < 180))
    {
        return angle;
    }

    Point newPoint = knee->getLateralPlateau() + tibiaImplantVector;

    if (sagital.eval(newPoint) > 0)
    {
        return angle;
    }
    else
    {
        return -angle;
    }
}

std::pair<Point, Point> ImplantsMatchFinalInfo::GetFemurBoneTEALine() const
{
    Line tLine(knee->getFemurVectorTEA(), knee->getLateralEpicondyle());
    Point medial = tLine.getProjectPoint(knee->getMedialEpicondyle());
    return std::make_pair(knee->getLateralEpicondyle(), medial);
}

std::pair<Point, Point> ImplantsMatchFinalInfo::GetFemurImplantTEALine() const
{
    double distance = (ImplantTools::getDistanceBetweenPoints(femurImplant.getPointP3(), femurImplant.getPointP4())) / 2.0;

    cv::Mat vectorimplantAPMat = femurRotation * (femurImplant.getDirectVectorTEA().ToMatPoint());
    Point vectorimplantAP = Point(vectorimplantAPMat);
    vectorimplantAP.normalice();

    cv::Mat pointImplantAPMat = (femurRotation * (femurImplant.getMidPlaneInterceptionPoint().ToMatPoint())) + femurTranslation;
    Point pointImplantAP = Point(pointImplantAPMat);

    Point a = pointImplantAP + distance * vectorimplantAP;
    Point b = pointImplantAP - distance * vectorimplantAP;

    return std::make_pair(a, b);
}

const itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::setFemurTEAAngle(double angle)
{
    double myAngle = (angle - GetFemurImplantTEAAngle()) * PI / 180.0;

    Point axisRotation = knee->getDirectVectorFemurAxis();

    if (knee->getIsRight() == true)
    {
        axisRotation = (-1.0) * axisRotation;
    }

    cv::Mat rotation = ImplantTools::getRotateMatrix(axisRotation, myAngle);
    femurRotation = rotation * femurRotation;

    return getITKFemurTransform();
}

double ImplantsMatchFinalInfo::GetFemurImplantTEAAngle() const
{
    Plane sagital;

    sagital.init(knee->getFemurVectorTEA(), knee->getFemurKneeCenter());

    if (sagital.eval(knee->getLateralEpicondyle()) < 0)
    {
        sagital.reverse();
    }

    Plane axial;
    axial.init(knee->getDirectVectorFemurAxis(), knee->getFemurKneeCenter());

    Point vectorBoneAP = knee->getFemurDirectVectorAP();
    cv::Mat vectorimplantAPMat = femurRotation * (femurImplant.getDirectVectorAP().ToMatPoint());
    Point vectorimplantAP = Point(vectorimplantAPMat);
    LegAngle legAngle;
    double angle = legAngle.getAngleBetweenVectors(axial.getProjectionVector(vectorBoneAP), axial.getProjectionVector(vectorimplantAP));

    if (!(angle > 0 && angle < 180))
    {
        return angle;
    }

    Point newPoint = knee->getFemurKneeCenter() + axial.getProjectionVector(vectorimplantAP);

    if (sagital.eval(newPoint) > 0)
    {
        return angle;
    }
    else
    {
        return -angle;
    }
}

const itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::setFemurPCAAngle(double angle)
{
    double myAngle = (angle - GetFemurImplantPCAAngle()) * PI / 180.0;

    Point axisRotation = knee->getDirectVectorFemurAxis();

    if (knee->getIsRight() == true)
    {
        axisRotation = (-1.0) * axisRotation;
    }

    cv::Mat rotation = ImplantTools::getRotateMatrix(axisRotation, myAngle);
    femurRotation = rotation * femurRotation;

    return getITKFemurTransform();
}

double ImplantsMatchFinalInfo::GetFemurImplantPCAAngle() const
{
    Plane axial;
    axial.init(knee->getDirectVectorFemurAxis(), knee->getFemurKneeCenter());

    if (axial.eval(knee->getHipCenter()) < 0)
    {
        axial.reverse();
    }

    Point vectorBonePCA = knee->getMedialCondyle() - knee->getLateralCondyle();
    vectorBonePCA = axial.getProjectionVector(vectorBonePCA);
    LegAngle legAngle;

    double sameDirection = 1.0;
    if (knee->getIsRight() == false)
    {
        sameDirection = -1.0;
    }

    cv::Mat vectorimplantTEAMat = femurRotation * (femurImplant.getDirectVectorTEA().ToMatPoint());
    Point vectorimplantTEA = Point(vectorimplantTEAMat);
    vectorimplantTEA = axial.getProjectionVector(vectorimplantTEA);
    vectorimplantTEA = sameDirection * vectorimplantTEA;
    
    double angle = legAngle.getAngleBetweenVectors(vectorBonePCA, vectorimplantTEA);

    if (!(angle > 0 && angle < 180))
    {
        return angle;
    }

    Plane femurHelp = axial.getPerpendicularPlane(knee->getMedialCondyle(), knee->getLateralCondyle());
    if (femurHelp.eval(knee->getFemurKneeCenter()) < 0)
    {
        femurHelp.reverse();
    }

    Point newPoint = knee->getLateralCondyle() + vectorimplantTEA;

    if (femurHelp.eval(newPoint) >= 0)
    {
        return angle;
    }
    else
    {
        return -angle;
    }
}


const itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::setFemurVarusAngle(double angle)
{
    double myAngle = (angle - GetFemurVarusAngle()) * PI / 180.0;

    Point axisRotation = knee->getFemurDirectVectorAP();
    if (knee->getIsRight() == false)
    {
        axisRotation = (-1.0) * axisRotation;
    }

    cv::Mat rotation = ImplantTools::getRotateMatrix(axisRotation, myAngle);
    femurRotation = rotation * femurRotation;

    return getITKFemurTransform();
}

double ImplantsMatchFinalInfo::GetFemurVarusAngle() const
{
    Point vectorBone = knee->getLateralEpicondyle() - knee->getMedialEpicondylePerp();
    LegAngle legAngle;
    double vectorSign = 1.0;

    if (knee->getIsRight() == true)
    {
        vectorSign = -1.0;
    }
    
    cv::Mat vectorImplantMat = femurRotation * (femurImplant.getDirectVectorTEA().ToMatPoint());
    Point vectorImplant = Point(vectorImplantMat);
    vectorImplant = vectorSign * vectorImplant;

    Plane coronal;
    coronal.init(knee->getFemurDirectVectorAP(), knee->getFemurKneeCenter());

    vectorBone = coronal.getProjectionVector(vectorBone);
    vectorImplant = coronal.getProjectionVector(vectorImplant);
    double angle = legAngle.getAngleBetweenVectors(vectorBone, vectorImplant);

    if (!(angle > 0 && angle < 180))
    {
        return angle;
    }

    Plane axial;
    axial.init(knee->getDirectVectorFemurAxis(), knee->getFemurKneeCenter());
    if (axial.eval(knee->getHipCenter()) > 0)
    {
        axial.reverse();
    }

    Point newPoint = knee->getFemurKneeCenter() + vectorImplant;
    if (axial.eval(newPoint) >= 0)
    {
        return angle;
    }
    else
    {
        return -angle;
    }
}

double ImplantsMatchFinalInfo::GetFemurFlexionAngle() const
{
    Plane sagital;
    sagital.init(knee->getFemurVectorTEA(), knee->getFemurKneeCenter());

    Plane axial;
    axial.init(knee->getDirectVectorFemurAxis(), knee->getFemurKneeCenter());
    if (axial.eval(knee->getHipCenter()) < 0)
    {
        axial.reverse();
    }

    Point boneAP = sagital.getProjectionVector(knee->getFemurDirectVectorAP());
    cv::Mat implantAPMat = femurRotation * (femurImplant.getDirectVectorAP().ToMatPoint());
    Point implantAP = Point(implantAPMat);
    implantAP = sagital.getProjectionVector(implantAP);
    LegAngle legAngle;
    double angle = legAngle.getAngleBetweenVectors(boneAP, implantAP);

    if (!(angle > 0 && angle < 180))
    {
        return angle;
    }

    Point newPoint = knee->getFemurKneeCenter() + implantAP;
    if (axial.eval(newPoint) >= 0)
    {
        return angle;
    }
    else
    {
        return -angle;
    }
}

const itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::setFemurFlexionAngle(double angle)
{
    double myAngle = (angle - GetFemurFlexionAngle()) * PI / 180.0;
    
    Point axisRotation;
    if (knee->getIsRight() == true)
    {
        axisRotation = knee->getLateralEpicondyle() - knee->getMedialEpicondylePerp();
    }
    else
    {
        axisRotation = knee->getMedialEpicondylePerp() - knee->getLateralEpicondyle();
    }

    cv::Mat rotation = ImplantTools::getRotateMatrix(axisRotation, myAngle);

    femurRotation = rotation * femurRotation;
    return getITKFemurTransform();
}

const itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::setTibiaVarusAngle(double angle)
{
    double myAngle = (angle - GetTibiaVarusAngle()) * PI / 180.0;

    Point axisRotation = knee->getTibiaDirectVectorAP();

    if (knee->getIsRight() == true)
    {
        axisRotation = (-1.0) * axisRotation;
    }

    cv::Mat rotation = ImplantTools::getRotateMatrix(axisRotation, myAngle);
    tibiaRotation = rotation * tibiaRotation;

    return getITKTibiaTransform();
}

double ImplantsMatchFinalInfo::GetTibiaVarusAngle() const
{
    Point crossVector = knee->getNormalVectorTibiaPlane().cross(knee->getTibiaDirectVectorAP());

    Plane sagital;
    sagital.init(crossVector, knee->getMedialPlateau());

    Point boneTEA = knee->getLateralPlateau() - sagital.getProjectionPoint(knee->getLateralPlateau());
    LegAngle legAngle;
    double vectorSign = 1.0;

    if (knee->getIsRight() == true)
    {
        vectorSign = -1.0;
    }

    cv::Mat implantTEAMat = tibiaRotation * (tibiaImplant.getTibiaVectorTEA().ToMatPoint());
    Point implantTEA = Point(implantTEAMat);

    implantTEA = vectorSign * implantTEA;

    Plane coronal;
    coronal.init(knee->getTibiaDirectVectorAP(), knee->getTibiaKneeCenter());

    boneTEA = coronal.getProjectionVector(boneTEA);
    implantTEA = coronal.getProjectionVector(implantTEA);

    
    double angle = legAngle.getAngleBetweenVectors(boneTEA, implantTEA);

    if (!(angle > 0 && angle < 180))
    {
        return angle;
    }

    Plane axial;
    axial.init(knee->getNormalVectorTibiaPlane(), knee->getTibiaKneeCenter());
    if (axial.eval(knee->getAnkleCenter()) > 0)
    {
        axial.reverse();
    }
    
    Point newPoint = knee->getTibiaKneeCenter() + implantTEA;

    if (axial.eval(newPoint) >= 0)
    {
        return angle;
    }
    else
    {
        return -angle;
    }
}

ResectionThickness ImplantsMatchFinalInfo::GetTibiaResection() const
{
    Plane planeTemp = tibiaImplant.getTibiaPlane();
    planeTemp.reverseByPoint(tibiaImplant.getExteriorPoint());

    Plane tibia = TransformPlane(planeTemp, tibiaRotation, tibiaTranslation);

    ResectionThickness result;
    result.lateral = tibia.eval(knee->getLateralPlateau());
    result.medial = tibia.eval(knee->getMedialPlateau());
    return result;
}

ResectionThickness ImplantsMatchFinalInfo::GetFemurResectionCoronal() const
{
    Plane planeTemp = femurImplant.getPlaneA();
    planeTemp.reverseByPoint(femurImplant.getCortexPoint(), false);

    Plane femur = TransformPlane(planeTemp, femurRotation, femurTranslation);

    ResectionThickness result;
    result.lateral = femur.eval(knee->getLateralCondyle());
    result.medial = femur.eval(knee->getMedialCondyle());
    return result;
}


itk::Vector< double, 3 > ImplantsMatchFinalInfo::SetThicknessTibiaMedial(double medialThickness)
{
    double finalThickness = medialThickness;// -tibiaCartilage;
    Point movePlateau = knee->getMedialPlateau();

    Plane planeTemp = tibiaImplant.getTibiaPlane();
    planeTemp.reverseByPoint(tibiaImplant.getExteriorPoint(), false);

    Plane tibia = TransformPlane(planeTemp, tibiaRotation, tibiaTranslation);
    Point oldPoint = tibia.getProjectionPoint(movePlateau);

    tibia.movePlane(movePlateau);
    tibia.movePlaneOnNormal(finalThickness);

    Point newMovePlateau = tibia.getProjectionPoint(movePlateau);

    tibiaTranslation = tibiaTranslation + (newMovePlateau.ToMatPoint() - oldPoint.ToMatPoint());

    Point newTranslationCV = Point(tibiaTranslation);

    itk::Vector< double, 3 > newTranslation;
    newTranslation[0] = newTranslationCV.x;
    newTranslation[1] = newTranslationCV.y;
    newTranslation[2] = newTranslationCV.z;

    return newTranslation;
}

itk::Vector< double, 3 > ImplantsMatchFinalInfo::SetThicknessTibiaLateral(double lateralThickness)
{
    double finalThickness = lateralThickness;// -tibiaCartilage;
    Point movePlateau = knee->getLateralPlateau();

    Plane planeTemp = tibiaImplant.getTibiaPlane();
    planeTemp.reverseByPoint(tibiaImplant.getExteriorPoint(), false);

    Plane tibia = TransformPlane(planeTemp, tibiaRotation, tibiaTranslation);
    Point oldPoint = tibia.getProjectionPoint(movePlateau);

    tibia.movePlane(movePlateau); 
    tibia.movePlaneOnNormal(finalThickness);

    Point newMovePlateau = tibia.getProjectionPoint(movePlateau);

    tibiaTranslation = tibiaTranslation + (newMovePlateau.ToMatPoint() - oldPoint.ToMatPoint());

    Point newTranslationCV = Point(tibiaTranslation);

    itk::Vector< double, 3 > newTranslation;
    newTranslation[0] = newTranslationCV.x;
    newTranslation[1] = newTranslationCV.y;
    newTranslation[2] = newTranslationCV.z;

    return newTranslation;
}

itk::Vector< double, 3 > ImplantsMatchFinalInfo::SetThicknessFemurCoronalMedial(double medialThickness)
{
    double finalThickness = medialThickness;// -femurCartilage;
    Point movePoint = knee->getMedialCondyle();

    Plane planeTemp = femurImplant.getPlaneA();
    planeTemp.reverseByPoint(femurImplant.getCortexPoint());

    Plane femur = TransformPlane(planeTemp, femurRotation, femurTranslation);

    Point oldPoint = femur.getProjectionPoint(movePoint);

    femur.movePlane(movePoint);
    femur.movePlaneOnNormal(finalThickness);

    Point newPoint = femur.getProjectionPoint(movePoint);

    femurTranslation = femurTranslation + (newPoint.ToMatPoint() - oldPoint.ToMatPoint());

    Point newTranslationCV = Point(femurTranslation);

    itk::Vector< double, 3 > newTranslation;
    newTranslation[0] = newTranslationCV.x;
    newTranslation[1] = newTranslationCV.y;
    newTranslation[2] = newTranslationCV.z;

    return newTranslation;
}

itk::Vector< double, 3 > ImplantsMatchFinalInfo::SetThicknessFemurCoronalLateral(double lateralThickness)
{
    double finalThickness = lateralThickness;// -femurCartilage;

    Point movePoint = knee->getLateralCondyle();

    Plane planeTemp = femurImplant.getPlaneA();
    planeTemp.reverseByPoint(femurImplant.getCortexPoint());

    Plane femur = TransformPlane(planeTemp, femurRotation, femurTranslation);

    Point oldPoint = femur.getProjectionPoint(movePoint);

    femur.movePlane(movePoint);
    femur.movePlaneOnNormal(finalThickness);

    Point newPoint = femur.getProjectionPoint(movePoint);

    femurTranslation = femurTranslation + (newPoint.ToMatPoint() - oldPoint.ToMatPoint());

    Point newTranslationCV = Point(femurTranslation);

    itk::Vector< double, 3 > newTranslation;
    newTranslation[0] = newTranslationCV.x;
    newTranslation[1] = newTranslationCV.y;
    newTranslation[2] = newTranslationCV.z;

    return newTranslation;
}

itk::Vector< double, 3 > ImplantsMatchFinalInfo::SetThicknessFemurAxialMedial(double medialThickness)
{
    double finalThickness = medialThickness;// -femurCartilage;
    Point movePoint = knee->getMedialInferiorFemurPoint();

    Plane planeTemp = femurImplant.getPlaneC();
    planeTemp.reverseByPoint(femurImplant.getCortexPoint());

    Plane femur = TransformPlane(planeTemp, femurRotation, femurTranslation);

    Point oldPoint = femur.getProjectionPoint(movePoint);

    femur.movePlane(movePoint);
    femur.movePlaneOnNormal(finalThickness);

    Point newPoint = femur.getProjectionPoint(movePoint);

    femurTranslation = femurTranslation + (newPoint.ToMatPoint() - oldPoint.ToMatPoint());

    Point newTranslationCV = Point(femurTranslation);

    itk::Vector< double, 3 > newTranslation;
    newTranslation[0] = newTranslationCV.x;
    newTranslation[1] = newTranslationCV.y;
    newTranslation[2] = newTranslationCV.z;

    return newTranslation;
}

itk::Vector< double, 3 > ImplantsMatchFinalInfo::SetThicknessFemurAxialLateral(double lateralThickness)
{
    double finalThickness = lateralThickness;// -femurCartilage;
    Point movePoint = knee->getLateralInferiorFemurPoint();

    Plane planeTemp = femurImplant.getPlaneC();
    planeTemp.reverseByPoint(femurImplant.getCortexPoint());

    Plane femur = TransformPlane(planeTemp, femurRotation, femurTranslation);

    Point oldPoint = femur.getProjectionPoint(movePoint);

    femur.movePlane(movePoint);
    femur.movePlaneOnNormal(finalThickness);

    Point newPoint = femur.getProjectionPoint(movePoint);

    femurTranslation = femurTranslation + (newPoint.ToMatPoint() - oldPoint.ToMatPoint());

    Point newTranslationCV = Point(femurTranslation);

    itk::Vector< double, 3 > newTranslation;
    newTranslation[0] = newTranslationCV.x;
    newTranslation[1] = newTranslationCV.y;
    newTranslation[2] = newTranslationCV.z;

    return newTranslation;
}

ResectionThickness ImplantsMatchFinalInfo::GetFemurResectionAxial() const
{
    Plane planeTemp = femurImplant.getPlaneC();
    planeTemp.reverseByPoint(femurImplant.getCortexPoint(), false);
    Plane femur = TransformPlane(planeTemp, femurRotation, femurTranslation);

    ResectionThickness result;
    result.lateral = femur.eval(knee->getLateralInferiorFemurPoint());
    result.medial = femur.eval(knee->getMedialInferiorFemurPoint());
    return result;
}

itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::getITKFemurTransform() const
{
    itk::Matrix< double, 3, 3 > rotation;
    itk::Vector< double, 3 > translate;

    rotation[0][0] = femurRotation.at<double>(0, 0);
    rotation[0][1] = femurRotation.at<double>(0, 1);
    rotation[0][2] = femurRotation.at<double>(0, 2);

    rotation[1][0] = femurRotation.at<double>(1, 0);
    rotation[1][1] = femurRotation.at<double>(1, 1);
    rotation[1][2] = femurRotation.at<double>(1, 2);

    rotation[2][0] = femurRotation.at<double>(2, 0);
    rotation[2][1] = femurRotation.at<double>(2, 1);
    rotation[2][2] = femurRotation.at<double>(2, 2);

    translate[0] = femurTranslation.at<double>(0, 0);
    translate[1] = femurTranslation.at<double>(1, 0);
    translate[2] = femurTranslation.at<double>(2, 0);

    itk::Rigid3DTransform<double>::Pointer transform = itk::VersorRigid3DTransform<double>::New();

    transform->SetMatrix(rotation);
    transform->SetOffset(translate);

    return transform;
}

itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::getITKTibiaTransform() const
{
    itk::Matrix< double, 3, 3 > rotation;
    itk::Vector< double, 3 > translate;

    rotation[0][0] = tibiaRotation.at<double>(0, 0);
    rotation[0][1] = tibiaRotation.at<double>(0, 1);
    rotation[0][2] = tibiaRotation.at<double>(0, 2);

    rotation[1][0] = tibiaRotation.at<double>(1, 0);
    rotation[1][1] = tibiaRotation.at<double>(1, 1);
    rotation[1][2] = tibiaRotation.at<double>(1, 2);

    rotation[2][0] = tibiaRotation.at<double>(2, 0);
    rotation[2][1] = tibiaRotation.at<double>(2, 1);
    rotation[2][2] = tibiaRotation.at<double>(2, 2);

    translate[0] = tibiaTranslation.at<double>(0, 0);
    translate[1] = tibiaTranslation.at<double>(1, 0);
    translate[2] = tibiaTranslation.at<double>(2, 0);

    itk::Rigid3DTransform<double>::Pointer transform = itk::VersorRigid3DTransform<double>::New();
    transform->SetMatrix(rotation);
    transform->SetOffset(translate);

    return transform;
}

cv::Mat ImplantsMatchFinalInfo::Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer transform) const
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

cv::Mat ImplantsMatchFinalInfo::Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer transform) const
{
    itk::Vector< double, 3 > translate = transform->GetOffset();
    cv::Mat result(3, 1, CV_64FC1);

    result.at <double>(0, 0) = translate[0];
    result.at <double>(1, 0) = translate[1];
    result.at <double>(2, 0) = translate[2];
    return result;
}

Plane ImplantsMatchFinalInfo::TransformPlane(const Plane& plane, const cv::Mat& rotation, const cv::Mat& translation) const
{
    cv::Mat transformNormalVector = rotation * plane.getNormalVectorMat();
    cv::Mat transformPoint = (rotation * plane.getPointMat()) + translation;
    Plane transformPlane;
    transformPlane.init(Point(transformNormalVector), Point(transformPoint));
    return transformPlane;
}

/*
itk::Vector< double, 3 > ImplantsMatchFinalInfo::SetThicknessFemurImplant(double lateralThickness, double medialThickness, bool alignByCortex)
{
    itk::Vector< double, 3 > newTranslation;

    Point oldTranslation = Point(femurTranslation);
    newTranslation[0] = oldTranslation.x;
    newTranslation[1] = oldTranslation.y;
    newTranslation[2] = oldTranslation.z;

    double newThickness = 0;
    if (knee->getIsVarus() == true)
    {
        newThickness = lateralThickness - femurLatThickness;
    }
    else
    {
        newThickness = medialThickness - femurMedialThickness;
    }

    Point distalA = knee->getFemurKneeCenter() + 1000.0 * knee->getFemurDirectVectorAP();
    Point distalC = knee->getHipCenter();

    Plane femurA = TransformPlane(femurImplant.getPlaneA(), femurRotation, femurTranslation);
    femurA.normalizeNormalVector();
    femurA.reverseByPoint(distalA);
    femurA.movePlaneOnNormal(newThickness);

    Plane femurMid = TransformPlane(femurImplant.getMidPlane(), femurRotation, femurTranslation);
    femurMid.normalizeNormalVector();

    Plane femurC = TransformPlane(femurImplant.getPlaneC(), femurRotation, femurTranslation);
    femurC.normalizeNormalVector();
    femurC.reverseByPoint(distalC);
    femurC.movePlaneOnNormal(newThickness);

    Plane femurE = TransformPlane(femurImplant.getPlaneE(), femurRotation, femurTranslation);
    femurE.normalizeNormalVector();

    cv::Mat pSeudoExpKneePointA = femurRotation * femurImplant.getPlaneA().getPointMat();
    cv::Mat pSeudoExpKneeMidPoint = femurRotation * femurImplant.getMidPlane().getPointMat();
    cv::Mat pSeudoExpKneePointC = femurRotation * femurImplant.getPlaneC().getPointMat();
    cv::Mat pSeudoExpKneePointE = femurRotation * femurImplant.getPlaneE().getPointMat();

    Point pSeudoKneePointA(pSeudoExpKneePointA);
    Point pSeudoKneeMidPoint(pSeudoExpKneeMidPoint);
    Point pSeudoKneePointC(pSeudoExpKneePointC);
    Point pSeudoKneePointE(pSeudoExpKneePointE);

    cv::Mat translationMatrix;
    bool result;

    if (alignByCortex == false)
    {
        std::vector<cv::Point3d> normalVectors;
        std::vector<double> biasVector;
        double bias = 0.0;

        normalVectors.push_back(femurA.getNormalVector().ToCVPoint());
        normalVectors.push_back(femurMid.getNormalVector().ToCVPoint());
        normalVectors.push_back(femurC.getNormalVector().ToCVPoint());

        cv::Mat A(normalVectors.size(), 3, CV_64F, normalVectors.data());

        bias = -(femurA.getBias() + (pSeudoKneePointA.dot(femurA.getNormalVector())));
        biasVector.push_back(bias);
        bias = -(femurMid.getBias() + (pSeudoKneeMidPoint.dot(femurMid.getNormalVector())));
        biasVector.push_back(bias);
        bias = -(femurC.getBias() + (pSeudoKneePointC.dot(femurC.getNormalVector())));
        biasVector.push_back(bias);
        cv::Mat B(biasVector.size(), 1, CV_64F, biasVector.data());
        result = cv::solve(A, B, translationMatrix);
    }
    else
    {
        std::vector<cv::Point3d> normalVectors;
        std::vector<double> biasVector;
        double bias = 0.0;

        normalVectors.push_back(femurE.getNormalVector().ToCVPoint());
        normalVectors.push_back(femurMid.getNormalVector().ToCVPoint());
        normalVectors.push_back(femurC.getNormalVector().ToCVPoint());
        cv::Mat A(normalVectors.size(), 3, CV_64F, normalVectors.data());
        bias = -(femurE.getBias() + (pSeudoKneePointE.dot(femurE.getNormalVector())));
        biasVector.push_back(bias);
        bias = -(femurMid.getBias() + (pSeudoKneeMidPoint.dot(femurMid.getNormalVector())));
        biasVector.push_back(bias);
        bias = -(femurC.getBias() + (pSeudoKneePointC.dot(femurC.getNormalVector())));
        biasVector.push_back(bias);
        cv::Mat B(biasVector.size(), 1, CV_64F, biasVector.data());
        result = cv::solve(A, B, translationMatrix);
    }

    if (result == true)
    {
        Point newTranslationMat = Point(translationMatrix);
        newTranslation[0] = newTranslationMat.x;
        newTranslation[1] = newTranslationMat.y;
        newTranslation[2] = newTranslationMat.z;
        femurTranslation = translationMatrix;
        femurLatThickness = lateralThickness;
        femurMedialThickness = medialThickness;
    }
    //std::cout << "Translation: " << translationMatrix << std::endl;
    return newTranslation;
}
*/

/*
itk::Vector< double, 3 > ImplantsMatchFinalInfo::SetThicknessTibiaImplant(double lateralThickness, double medialThickness)
{
    double newThickness = 0;
    Point movePlateau;
    if (knee->getIsVarus() == true)
    {
        newThickness = lateralThickness;
        movePlateau = knee->getLateralPlateau();
    }
    else
    {
        newThickness = medialThickness;
        movePlateau = knee->getMedialPlateau();
    }

    Point distalPoint = knee->getAnkleCenter();

    Plane tibia = TransformPlane(tibiaImplant.getTibiaPlane(), tibiaRotation, tibiaTranslation);
    Point oldPoint = tibia.getProjectionPoint(movePlateau);
    tibia.movePlane(movePlateau);
    tibia.reverseByPoint(distalPoint);
    tibia.movePlaneOnNormal(newThickness);

    Point newMovePlateau = tibia.getProjectionPoint(movePlateau);

    tibiaTranslation = tibiaTranslation + (newMovePlateau.ToMatPoint() - oldPoint.ToMatPoint());

    Point newTranslationCV = Point(tibiaTranslation);

    itk::Vector< double, 3 > newTranslation;
    newTranslation[0] = newTranslationCV.x;
    newTranslation[1] = newTranslationCV.y;
    newTranslation[2] = newTranslationCV.z;

    return newTranslation;
}
*/


void ImplantsMatchFinalInfo::SetCartilageFemur(double pCartilage)
{
    femurCartilage = pCartilage;
    knee->setFemurCartilage(pCartilage);
}

void ImplantsMatchFinalInfo::SetCartilageTibia(double pCartilage)
{
    tibiaCartilage = pCartilage;
    knee->setTibiaCartilage(pCartilage);
}


void ImplantsMatchFinalInfo::setFemurTransform(const itk::Rigid3DTransform<>::Pointer pImplantToBoneFemurTransform)
{
    femurRotation = Rigid3DTransformToCVRotation(pImplantToBoneFemurTransform);
    femurTranslation = Rigid3DTransformToCVTranslation(pImplantToBoneFemurTransform);
    updateFemurImplantVectors();
}

void ImplantsMatchFinalInfo::setTibiaTransform(const itk::Rigid3DTransform<>::Pointer pImplantToBoneTibiaTransform)
{
    tibiaRotation = Rigid3DTransformToCVRotation(pImplantToBoneTibiaTransform);
    tibiaTranslation = Rigid3DTransformToCVTranslation(pImplantToBoneTibiaTransform);
    updateTibiaImplantVectors();
}

void ImplantsMatchFinalInfo::test()
{
    ResectionThickness tibia = GetTibiaResection();
    ResectionThickness axialFemur = GetFemurResectionAxial();
    ResectionThickness coronalFemur = GetFemurResectionCoronal();
    double varusAngle = GetFemurVarusAngle();
    double flexionAngle = GetFemurFlexionAngle();
    double slopeAngle = GetTibiaImplantSlopeAngle();
    double tibiaRotationAngle = GetTibiaImplantRotationAngle();
    double TEAAngle = GetFemurImplantTEAAngle();
    double PCAAngle = GetFemurImplantPCAAngle();
    double tibiaVarus = GetTibiaVarusAngle();

    std::cout << "Femur Flexion: " << flexionAngle << std::endl;
    std::cout << "Femur Varus: " << varusAngle << std::endl;
    std::cout << "Tibia Slope: " << slopeAngle << std::endl;
    std::cout << "Tibia rotation: " << tibiaRotationAngle << std::endl;
    std::cout << "Tibia Varus: " << tibiaVarus << std::endl;
    std::cout << "TEA: " << TEAAngle << std::endl;
    std::cout << "PCA: " << PCAAngle << std::endl;
    std::cout << "Tibia resection lateral: " << tibia.lateral << " Medial: " << tibia.medial << std::endl;
    std::cout << "Femur resection axial lateral: " << axialFemur.lateral << " Medial: " << axialFemur.medial << std::endl;
    std::cout << "Femur resection coronal lateral: " << coronalFemur.lateral << " Medial: " << coronalFemur.medial << std::endl;

    std::cout << "Translation femur: " << femurTranslation << std::endl;
    std::cout << "Translation tibia: " << tibiaTranslation << std::endl;

    Point onBone = GetFemurBoneTEALine().first - GetFemurBoneTEALine().second;
    Point onImplant = GetFemurImplantTEALine().first - GetFemurImplantTEALine().second;
    onBone.normalice();
    onImplant.normalice();

    std::cout << "femur bone tea vector: " << onBone << std::endl;
    std::cout << "femur implant tea vector: " << onImplant << std::endl;
}