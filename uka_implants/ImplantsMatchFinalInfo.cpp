#include "ImplantsMatchFinalInfo.hpp"
#include "LegAngle.hpp"
#include "ImplantTools.hpp"
#include <itkVersorRigid3DTransform.h>

using namespace UKA::IMPLANTS;

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
    //this->boneKneeCapPath = pKnee->getKneeGroovePath();
    //this->implantKneeCapPath = pFemurImplant.GetKneeCapPath();

    updateFemurImplantVectors();
    updateTibiaImplantVectors();
}

void ImplantsMatchFinalInfo::updateFemurImplantVectors()
{
    cv::Mat axisTemp = femurRotation * (femurImplant.getDirectVectorFemurAxis().ToMatPoint());
    cv::Mat apTemp = femurRotation * (femurImplant.getDirectVectorAP().ToMatPoint());
    cv::Mat TEATemp = femurRotation * (femurImplant.getDirectVectorTEA().ToMatPoint());
    cv::Mat middleTemp = femurRotation * (femurImplant.getRodTopPointProjectedOnBase().ToMatPoint()) + femurTranslation;

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

itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::FemurImplantToTibiaImplant()
{
	std::vector<cv::Point3d> femurVectors;
	std::vector<cv::Point3d> tibiaVectors;

	Point implantFemurAxis = femurImplant.getDirectVectorFemurAxis();
	Point implantFemurAP = femurImplant.getDirectVectorAP();
	Point implantFemurTEA = implantFemurAxis.cross(implantFemurAP);
	implantFemurTEA.normalice();

	femurVectors.push_back(implantFemurTEA.ToCVPoint());
	femurVectors.push_back(implantFemurAxis.ToCVPoint());
	femurVectors.push_back(implantFemurAP.ToCVPoint());

	cv::Mat implantTibiaAxisMat = tibiaRotation * tibiaImplant.getTibiaNormalVector().ToMatPoint();
	cv::Mat implantTibiaAPMat = tibiaRotation * tibiaImplant.getTibiaVectorAP().ToMatPoint();

	Point implantTibiaAxis = Point(implantTibiaAxisMat);
	Point implantTibiaAP = Point(implantTibiaAPMat);
	Point implantTibiaTEA = implantTibiaAxis.cross(implantTibiaAP);
	implantTibiaTEA.normalice();

	tibiaVectors.push_back(implantTibiaTEA.ToCVPoint());
	tibiaVectors.push_back(implantTibiaAxis.ToCVPoint());
	tibiaVectors.push_back(implantTibiaAP.ToCVPoint());

	cv::Mat implantFemurMatrix = cv::Mat(femurVectors.size(), 3, CV_64F, femurVectors.data());
	cv::Mat implantTibiaMatrix = cv::Mat(tibiaVectors.size(), 3, CV_64F, tibiaVectors.data());

	cv::Mat inverse = (implantFemurMatrix.t()).inv();
	femurRotation = (implantTibiaMatrix.t()) * inverse;

	cv::Mat refPointMat = tibiaRotation * tibiaImplant.getPlateauRefPointUp().ToMatPoint() + tibiaTranslation;
	femurTranslation = refPointMat - femurRotation * femurImplant.getRodTopPointProjectedOnBaseExterior().ToMatPoint();

	return getITKFemurTransform();
}

itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::TibiaImplantToFemurImplant()
{
	std::vector<cv::Point3d> femurVectors;
	std::vector<cv::Point3d> tibiaVectors;

	Point implantTibiaAxis = tibiaImplant.getTibiaNormalVector();
	Point implantTibiaAP = tibiaImplant.getTibiaVectorAP();
	Point implantTibiaTEA = implantTibiaAxis.cross(implantTibiaAP);
	implantTibiaTEA.normalice();

	tibiaVectors.push_back(implantTibiaTEA.ToCVPoint());
	tibiaVectors.push_back(implantTibiaAxis.ToCVPoint());
	tibiaVectors.push_back(implantTibiaAP.ToCVPoint());

	cv::Mat implantFemurAxisMat = femurRotation * femurImplant.getDirectVectorFemurAxis().ToMatPoint();
	cv::Mat implantFemurAPMat = femurRotation * femurImplant.getDirectVectorAP().ToMatPoint();

	Point implantFemurAxis = Point(implantFemurAxisMat);
	Point implantFemurAP = Point(implantFemurAPMat);
	Point implantFemurTEA = implantFemurAxis.cross(implantFemurAP);
	implantFemurTEA.normalice();

	femurVectors.push_back(implantFemurTEA.ToCVPoint());
	femurVectors.push_back(implantFemurAxis.ToCVPoint());
	femurVectors.push_back(implantFemurAP.ToCVPoint());

	cv::Mat implantTibiaMatrix = cv::Mat(tibiaVectors.size(), 3, CV_64F, tibiaVectors.data());
	cv::Mat implantFemurMatrix = cv::Mat(femurVectors.size(), 3, CV_64F, femurVectors.data());

	cv::Mat inverse = (implantTibiaMatrix.t()).inv();
	tibiaRotation = (implantFemurMatrix.t()) * inverse;

	cv::Mat refPointMat = femurRotation * femurImplant.getRodTopPointProjectedOnBaseExterior().ToMatPoint() + femurTranslation;
	tibiaTranslation = refPointMat - tibiaRotation * tibiaImplant.getPlateauRefPointUp().ToMatPoint();

	return getITKTibiaTransform();
}

const itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::setFemurImplant(const FemurImplant& pFemurImplant)
{
    femurImplant = pFemurImplant;
    //this->implantKneeCapPath = pFemurImplant.GetKneeCapPath();

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
    femurTranslation = femurImplantMiddlePoint.ToMatPoint() - (femurRotation * femurImplant.getRodTopPointProjectedOnBase().ToMatPoint());

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

/*
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
        cv::Mat newPoint = femurRotation * (implantKneeCapPath[i].ToMatPoint()) + femurTranslation;
        Point proj = Point(newPoint);
        pointsOut.push_back(proj.ToITKPoint());
    }

    return pointsOut;
}
*/

const itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::setTibiaSlopeAngle(double angle)
{
    double myAngle = angle * PI / 180.0;

	Point tibiaAxis = knee->getTibiaKneeCenter() - knee->getAnkleCenter();
	Plane tibiaHelp;
	tibiaHelp.init(tibiaAxis, knee->getTibiaKneeCenter());

	Point boneAP = -(knee->getTibiaTubercle() - knee->getPclCenterPoint());
	boneAP = tibiaHelp.getProjectionVector(boneAP);
	boneAP.normalice();

    Point axisRotation = tibiaHelp.getNormalVector().cross(boneAP);

	////////////////////////////////////////////////////////// Rotation

	cv::Mat rotation = ImplantTools::getRotateMatrix(axisRotation, myAngle);

	auto mainVectorMat = rotation * boneAP.ToMatPoint();
	auto mainVector = Point(mainVectorMat);

	Plane projectionPlane;
	projectionPlane.init(axisRotation, knee->getTibiaKneeCenter());

	Plane implantPlane;
	Point axisRotationImplant = tibiaImplant.getTibiaVectorTEA();
	implantPlane.init(axisRotationImplant, tibiaImplant.getCentralPoint());
	implantPlane.transformPlane(tibiaRotation, tibiaTranslation);

	Point newVectorFromImplant = ImplantTools::getOriginalVectorFromProjectionWithPlanes(projectionPlane, mainVector, implantPlane);

	auto baseVectorFromImplantMat = tibiaRotation * tibiaImplant.getTibiaVectorAP().ToMatPoint();
	Point baseVectorFromImplant = -Point(baseVectorFromImplantMat);

	myAngle = ImplantTools::getAngleBetweenVectors(newVectorFromImplant, baseVectorFromImplant);

	rotation = ImplantTools::getRotateMatrix(baseVectorFromImplant.cross(newVectorFromImplant), myAngle);

	//////////////////////////////////////////////////////////////////////

	tibiaRotation = rotation * tibiaRotation;

    knee->setTibiaSlope(angle);

    return getITKTibiaTransform();
}

double ImplantsMatchFinalInfo::GetTibiaImplantSlopeAngle() const
{
    Point tibiaAxis = knee->getTibiaKneeCenter() - knee->getAnkleCenter();
    Plane tibiaHelp;
    tibiaHelp.init(tibiaAxis, knee->getTibiaKneeCenter());

    Point boneAP = -(knee->getTibiaTubercle() - knee->getPclCenterPoint());
    boneAP = tibiaHelp.getProjectionVector(boneAP);

    Point implantAP = -(tibiaImplant.getTibiaVectorAP());
    cv::Mat implantAPMat = tibiaRotation * (implantAP.ToMatPoint());
    implantAP = Point(implantAPMat);

    Plane sagital;
	sagital.init(tibiaAxis.cross(boneAP), knee->getTibiaKneeCenter());

    if (tibiaHelp.eval(knee->getAnkleCenter()) < 0)
    {
        tibiaHelp.reverse();
    }

    boneAP = sagital.getProjectionVector(boneAP);
    implantAP = sagital.getProjectionVector(implantAP);

    LegAngle legAngle;
    double angle = legAngle.getAngleBetweenVectors(boneAP, implantAP);

    Point newPoint = knee->getTibiaKneeCenter() + implantAP;

    if (tibiaHelp.eval(newPoint) >= 0)
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
    double myAngle = angle * PI / 180.0;

	Point tibiaAxis = knee->getTibiaKneeCenter() - knee->getAnkleCenter();
	Plane tibiaHelp;
	tibiaHelp.init(tibiaAxis, knee->getTibiaKneeCenter());

	Point boneAP = knee->getTibiaTubercle() - knee->getPclCenterPoint();;
	boneAP = tibiaHelp.getProjectionVector(boneAP);

    Point axisRotation = tibiaHelp.getNormalVector();

    if (knee->getIsRight() == true)
    {
        axisRotation = -axisRotation;
    }

	////////////////////////////////////////////////////////// Rotation

	cv::Mat rotation = ImplantTools::getRotateMatrix(axisRotation, myAngle);

	auto mainVectorMat = rotation * boneAP.ToMatPoint();
	auto mainVector = Point(mainVectorMat);

	Plane projectionPlane;
	projectionPlane.init(axisRotation, knee->getTibiaKneeCenter());

	Plane implantPlane;
	Point axisRotationImplant = tibiaImplant.getTibiaNormalVector();
	implantPlane.init(axisRotationImplant, tibiaImplant.getCentralPoint());
	implantPlane.transformPlane(tibiaRotation, tibiaTranslation);

	Point newVectorFromImplant = ImplantTools::getOriginalVectorFromProjectionWithPlanes(projectionPlane, mainVector, implantPlane);

	auto baseVectorFromImplantMat = tibiaRotation * tibiaImplant.getTibiaVectorAP().ToMatPoint();
	Point baseVectorFromImplant = Point(baseVectorFromImplantMat);

	myAngle = ImplantTools::getAngleBetweenVectors(newVectorFromImplant, baseVectorFromImplant);

	rotation = ImplantTools::getRotateMatrix(baseVectorFromImplant.cross(newVectorFromImplant), myAngle);

	//////////////////////////////////////////////////////////////////////

    tibiaRotation = rotation * tibiaRotation;

    return getITKTibiaTransform();
}

double ImplantsMatchFinalInfo::GetTibiaImplantRotationAngle() const
{
	Point tibiaAxis = knee->getTibiaKneeCenter() - knee->getAnkleCenter();
	Plane transversal;
	transversal.init(tibiaAxis, knee->getTibiaKneeCenter());

	Point boneAP = knee->getTibiaTubercle() - knee->getPclCenterPoint();
	boneAP = transversal.getProjectionVector(boneAP);
	boneAP.normalice();

	Point vectorLateralTEA = transversal.getNormalVector().cross(boneAP);

	if (knee->getIsRight() == true)
	{
		vectorLateralTEA = -vectorLateralTEA;
	}

    Plane sagital;
    sagital.init(vectorLateralTEA, knee->getTibiaKneeCenter());

    cv::Mat tibiaImplantVectorMat = tibiaRotation * (tibiaImplant.getTibiaVectorAP().ToMatPoint());
    Point tibiaImplantVector = Point(tibiaImplantVectorMat);
    tibiaImplantVector = transversal.getProjectionVector(tibiaImplantVector);
    LegAngle legAngle;

    double angle = legAngle.getAngleBetweenVectors(boneAP, tibiaImplantVector);

    Point newPoint = knee->getTibiaKneeCenter() + tibiaImplantVector;

    if (sagital.eval(newPoint) >= 0)
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
    double distance = (femurImplant.getWidthSize()) / 2.0;
	Point tCentral = femurImplant.getMidPlane().getProjectionPoint(femurImplant.getRodBasePoint());

    cv::Mat vectorimplantTEAMat = femurRotation * (femurImplant.getDirectVectorTEA().ToMatPoint());
    Point vectorimplantTEA = Point(vectorimplantTEAMat);
    vectorimplantTEA.normalice();

	cv::Mat tCentralTransform = femurRotation * (tCentral.ToMatPoint()) + femurTranslation;
	tCentral = Point(tCentralTransform);

    Point a = tCentral + distance * vectorimplantTEA;
    Point b = tCentral - distance * vectorimplantTEA;

    return std::make_pair(a, b);
}

const itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::setFemurTEAAngle(double angle)
{
    double myAngle = angle * PI / 180.0;

    Point axisRotation = knee->getDirectVectorFemurAxis();

    if (knee->getIsRight() == true)
    {
        axisRotation = -axisRotation;
    }

	////////////////////////////////////////////////////////// Rotation

	cv::Mat rotation = ImplantTools::getRotateMatrix(axisRotation, myAngle);

	auto mainVectorMat = rotation * knee->getFemurVectorTEA().ToMatPoint();
	auto mainVector = Point(mainVectorMat);

	Plane projectionPlane;
	projectionPlane.init(axisRotation, knee->getFemurKneeCenter());

	Plane implantPlane;
	Point axisRotationImplant = femurImplant.getDirectVectorFemurAxis();
	implantPlane.init(axisRotationImplant, femurImplant.getRodTopPointProjectedOnBase());
	implantPlane.transformPlane(femurRotation, femurTranslation);

	Point newVectorFromImplant = ImplantTools::getOriginalVectorFromProjectionWithPlanes(projectionPlane, mainVector, implantPlane);
	
	auto baseVectorFromImplantMat = femurRotation * femurImplant.getDirectVectorTEA().ToMatPoint();
	Point baseVectorFromImplant = Point(baseVectorFromImplantMat);

	myAngle = ImplantTools::getAngleBetweenVectors(newVectorFromImplant, baseVectorFromImplant);

	rotation = ImplantTools::getRotateMatrix(baseVectorFromImplant.cross(newVectorFromImplant), myAngle);

	//////////////////////////////////////////////////////////////////////

    femurRotation = rotation * femurRotation;

    return getITKFemurTransform();
}

double ImplantsMatchFinalInfo::GetFemurImplantTEAAngle() const
{
    Plane coronal;
	coronal.init(knee->getFemurDirectVectorAP(), knee->getFemurKneeCenter());

    Plane axial;
    axial.init(knee->getDirectVectorFemurAxis(), knee->getFemurKneeCenter());

    Point vectorBoneTEA = knee->getFemurVectorTEA();
    cv::Mat vectorimplantTEAMat = femurRotation * (femurImplant.getDirectVectorTEA().ToMatPoint());
    Point vectorimplantTEA = Point(vectorimplantTEAMat);
    LegAngle legAngle;
    double angle = legAngle.getAngleBetweenVectors(axial.getProjectionVector(vectorBoneTEA), axial.getProjectionVector(vectorimplantTEA));

	Point refVector = axial.getProjectionVector(vectorimplantTEA);

	if (knee->getIsRight() == false)
	{
		refVector = -refVector;
	}

    Point newPoint = knee->getFemurKneeCenter() + refVector;

    if (coronal.eval(newPoint) > 0)
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
	double myAngle = angle * PI / 180.0;

	Point axisRotation = knee->getDirectVectorFemurAxis();

	if (knee->getIsRight() == true)
	{
		axisRotation = -axisRotation;
	}

	Point ref = knee->getMedialCondyle() - knee->getLateralCondyle();
	Plane projectionPlane;
	projectionPlane.init(axisRotation, knee->getFemurKneeCenter());
	ref = projectionPlane.getProjectionVector(ref);
	ref.normalice();

	////////////////////////////////////////////////////////// Rotation

	cv::Mat rotation = ImplantTools::getRotateMatrix(axisRotation, myAngle);

	auto mainVectorMat = rotation * ref.ToMatPoint();
	auto mainVector = Point(mainVectorMat);

	Plane implantPlane;
	Point axisRotationImplant = femurImplant.getDirectVectorFemurAxis();
	implantPlane.init(axisRotationImplant, femurImplant.getRodTopPointProjectedOnBase());
	implantPlane.transformPlane(femurRotation, femurTranslation);

	Point newVectorFromImplant = ImplantTools::getOriginalVectorFromProjectionWithPlanes(projectionPlane, mainVector, implantPlane);

	auto baseVectorFromImplantMat = femurRotation * femurImplant.getDirectVectorTEA().ToMatPoint();
	Point baseVectorFromImplant = Point(baseVectorFromImplantMat);

	if (knee->getIsRight() == false)
	{
		baseVectorFromImplant = -baseVectorFromImplant;
	}

	myAngle = ImplantTools::getAngleBetweenVectors(newVectorFromImplant, baseVectorFromImplant);

	rotation = ImplantTools::getRotateMatrix(baseVectorFromImplant.cross(newVectorFromImplant), myAngle);

	//////////////////////////////////////////////////////////////////////

	femurRotation = rotation * femurRotation;

	/*auto myAP = femurRotation * femurImplant.getDirectVectorAP().ToMatPoint();
	auto myTEA = femurRotation * femurImplant.getDirectVectorTEA().ToMatPoint();
	
	if (knee->getIsRight() == false)
	{
		myTEA = -myTEA;
	}

	auto projAP = axial.getProjectionVector(Point(myAP));
	auto projTEA = axial.getProjectionVector(Point(myTEA));
	
	std::cout << "TEA: " << ImplantTools::getAngleBetweenVectorsDegree(projTEA, vectorTea) << std::endl;
	std::cout << "PCA: " << ImplantTools::getAngleBetweenVectorsDegree(projTEA, vectorBonePCA) << std::endl;
	std::cout << "PCA oficial: " << ImplantTools::getAngleBetweenVectorsDegree(vectorTea, vectorBonePCA) << std::endl;
	std::cout << "Implant AP and TEA 3: " << ImplantTools::getAngleBetweenVectorsDegree(projAP, projTEA) << std::endl;
	std::cout << "Implant AP and TEA 1: " << ImplantTools::getAngleBetweenVectorsDegree(femurImplant.getDirectVectorAP(), femurImplant.getDirectVectorTEA()) << std::endl;
	std::cout << "Implant AP and TEA 2: " << ImplantTools::getAngleBetweenVectorsDegree(Point(myAP), Point(myTEA)) << std::endl;*/

	return getITKFemurTransform();






    //Point axisRotation = knee->getDirectVectorFemurAxis();

    //if (knee->getIsRight() == true)
    //{
    //    axisRotation = -axisRotation;
    //}

    //cv::Mat rotation = ImplantTools::getRotateMatrix(axisRotation, myAngle);
    //femurRotation = rotation * femurRotation;

    //return getITKFemurTransform();

	//return setFemurTEAAngle(myAngle);
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
    double myAngle = angle * PI / 180.0;

    Point axisRotation = knee->getFemurDirectVectorAP();

    if (knee->getIsRight() == false)
    {
        axisRotation = -axisRotation;
    }

	////////////////////////////////////////////////////////// Rotation

	cv::Mat rotation = ImplantTools::getRotateMatrix(axisRotation, myAngle);

	auto mainVectorMat = rotation * knee->getDirectVectorFemurAxis().ToMatPoint();
	auto mainVector = Point(mainVectorMat);

	Plane projectionPlane;
	projectionPlane.init(axisRotation, knee->getFemurKneeCenter());

	Plane implantPlane;
	Point axisRotationImplant = femurImplant.getDirectVectorAP();
	implantPlane.init(axisRotationImplant, femurImplant.getRodTopPointProjectedOnBase());
	implantPlane.transformPlane(femurRotation, femurTranslation);

	Point newVectorFromImplant = ImplantTools::getOriginalVectorFromProjectionWithPlanes(projectionPlane, mainVector, implantPlane);

	auto baseVectorFromImplantMat = femurRotation * femurImplant.getDirectVectorFemurAxis().ToMatPoint();
	Point baseVectorFromImplant = Point(baseVectorFromImplantMat);

	myAngle = ImplantTools::getAngleBetweenVectors(newVectorFromImplant, baseVectorFromImplant);

	rotation = ImplantTools::getRotateMatrix(baseVectorFromImplant.cross(newVectorFromImplant), myAngle);

	//////////////////////////////////////////////////////////////////////

    femurRotation = rotation * femurRotation;

    return getITKFemurTransform();
}

double ImplantsMatchFinalInfo::GetFemurVarusAngle() const
{
    Point vectorBone = knee->getDirectVectorFemurAxis();
    LegAngle legAngle;
    
    cv::Mat vectorImplantMat = femurRotation * (femurImplant.getDirectVectorFemurAxis().ToMatPoint());
    Point vectorImplant = Point(vectorImplantMat);

    Plane coronal;
    coronal.init(knee->getFemurDirectVectorAP(), knee->getFemurKneeCenter());

    vectorImplant = coronal.getProjectionVector(vectorImplant);
    double angle = legAngle.getAngleBetweenVectors(vectorBone, vectorImplant);

    Plane sagital;
	sagital.init(knee->getFemurVectorLateralTEA(), knee->getFemurKneeCenter());

    Point newPoint = knee->getFemurKneeCenter() + vectorImplant;
    if (sagital.eval(newPoint) >= 0)
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

    Point boneAP = knee->getFemurDirectVectorAP();
    cv::Mat implantAPMat = femurRotation * (femurImplant.getDirectVectorAP().ToMatPoint());
    Point implantAP = Point(implantAPMat);
    implantAP = sagital.getProjectionVector(implantAP);
    LegAngle legAngle;
    double angle = legAngle.getAngleBetweenVectors(boneAP, implantAP);

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
	double myAngle = angle * PI / 180.0;

	Point axisRotation;
	if (knee->getIsRight() == true)
	{
		axisRotation = knee->getLateralEpicondyle() - knee->getMedialEpicondylePerp();
	}
	else
	{
		axisRotation = knee->getMedialEpicondylePerp() - knee->getLateralEpicondyle();
	}
	axisRotation.normalice();

	////////////////////////////////////////////////////////// Rotation

	cv::Mat rotation = ImplantTools::getRotateMatrix(axisRotation, myAngle);

	auto mainVectorMat = rotation * knee->getFemurDirectVectorAP().ToMatPoint();
	auto mainVector = Point(mainVectorMat);

	Plane projectionPlane;
	projectionPlane.init(axisRotation, knee->getFemurKneeCenter());

	Plane implantPlane;
	Point axisRotationImplant = femurImplant.getDirectVectorTEA();
	implantPlane.init(axisRotationImplant, femurImplant.getRodTopPointProjectedOnBase());
	implantPlane.transformPlane(femurRotation, femurTranslation);

	Point newVectorFromImplant = ImplantTools::getOriginalVectorFromProjectionWithPlanes(projectionPlane, mainVector, implantPlane);

	auto baseVectorFromImplantMat = femurRotation * femurImplant.getDirectVectorAP().ToMatPoint();
	Point baseVectorFromImplant = Point(baseVectorFromImplantMat);

	myAngle = ImplantTools::getAngleBetweenVectors(newVectorFromImplant, baseVectorFromImplant);

	rotation = ImplantTools::getRotateMatrix(baseVectorFromImplant.cross(newVectorFromImplant), myAngle);

	//////////////////////////////////////////////////////////////////////

    femurRotation = rotation * femurRotation;
    return getITKFemurTransform();
}

const itk::Rigid3DTransform<>::Pointer ImplantsMatchFinalInfo::setTibiaVarusAngle(double angle)
{
	double myAngle = angle * PI / 180.0;

	Point tibiaAxis = knee->getTibiaKneeCenter() - knee->getAnkleCenter();
	Plane tibiaHelp;
	tibiaHelp.init(tibiaAxis, knee->getTibiaKneeCenter());

	Point boneAP = knee->getTibiaTubercle() - knee->getPclCenterPoint();
	boneAP = tibiaHelp.getProjectionVector(boneAP);

	Point axisRotation = boneAP;

	if (knee->getIsRight() == true)
	{
		axisRotation = -axisRotation;
	}

	////////////////////////////////////////////////////////// Rotation

	cv::Mat rotation = ImplantTools::getRotateMatrix(axisRotation, myAngle);

	auto mainVectorMat = rotation * tibiaHelp.getNormalVector().ToMatPoint();
	auto mainVector = Point(mainVectorMat);

	Plane projectionPlane;
	projectionPlane.init(axisRotation, knee->getTibiaKneeCenter());

	Plane implantPlane;
	Point axisRotationImplant = tibiaImplant.getTibiaVectorAP();
	implantPlane.init(axisRotationImplant, tibiaImplant.getCentralPoint());
	implantPlane.transformPlane(tibiaRotation, tibiaTranslation);

	Point newVectorFromImplant = ImplantTools::getOriginalVectorFromProjectionWithPlanes(projectionPlane, mainVector, implantPlane);

	auto baseVectorFromImplantMat = tibiaRotation * tibiaImplant.getTibiaNormalVector().ToMatPoint();
	Point baseVectorFromImplant = Point(baseVectorFromImplantMat);

	myAngle = ImplantTools::getAngleBetweenVectors(newVectorFromImplant, baseVectorFromImplant);

	rotation = ImplantTools::getRotateMatrix(baseVectorFromImplant.cross(newVectorFromImplant), myAngle);

	//////////////////////////////////////////////////////////////////////

    tibiaRotation = rotation * tibiaRotation;

    return getITKTibiaTransform();
}

double ImplantsMatchFinalInfo::GetTibiaVarusAngle() const
{
	Point tibiaAxis = knee->getTibiaKneeCenter() - knee->getAnkleCenter();
	Plane tibiaHelp;
	tibiaHelp.init(tibiaAxis, knee->getTibiaKneeCenter());

	Point boneAP = knee->getTibiaTubercle() - knee->getPclCenterPoint();
	boneAP = tibiaHelp.getProjectionVector(boneAP);
	boneAP.normalice();

	Point vectorMedialTEA = tibiaHelp.getNormalVector().cross(boneAP);

	if (knee->getIsRight() == false)
	{
		vectorMedialTEA = -vectorMedialTEA;
	}

    Plane sagital;
    sagital.init(vectorMedialTEA, knee->getTibiaKneeCenter());

	Plane coronal;
	coronal.init(boneAP, knee->getTibiaKneeCenter());

	cv::Mat implantAxisMat = tibiaRotation * (tibiaImplant.getTibiaNormalVector().ToMatPoint());
	Point implantAxis = Point(implantAxisMat);
	Point boneAxis = tibiaHelp.getNormalVector();

	boneAxis = coronal.getProjectionVector(boneAxis);
	implantAxis = coronal.getProjectionVector(implantAxis);

    LegAngle legAngle;

    double angle = legAngle.getAngleBetweenVectors(boneAxis, implantAxis);
    
    Point newPoint = knee->getTibiaKneeCenter() + implantAxis;

    if (sagital.eval(newPoint) >= 0)
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

	if (knee->getSurgerySide() == SurgerySideEnum::KLateral)
	{
		result.SurgerySide = tibia.eval(knee->getLateralPlateau());
		result.NoSurgerySide = tibia.eval(knee->getMedialPlateau());
	}
	else
	{
		result.NoSurgerySide = tibia.eval(knee->getLateralPlateau());
		result.SurgerySide = tibia.eval(knee->getMedialPlateau());
	}

    return result;
}

ResectionThickness ImplantsMatchFinalInfo::GetFemurResectionCoronal() const
{
    Plane planeTemp = femurImplant.getPosterior();
    planeTemp.reverseByPoint(femurImplant.getRodBasePoint(), false);

    Plane femur = TransformPlane(planeTemp, femurRotation, femurTranslation);

    ResectionThickness result;

	if (knee->getSurgerySide() == SurgerySideEnum::KLateral)
	{
		result.SurgerySide = femur.eval(knee->getLateralCondyle());
		result.NoSurgerySide = femur.eval(knee->getMedialCondyle());
	}
	else
	{
		result.NoSurgerySide = femur.eval(knee->getLateralCondyle());
		result.SurgerySide = femur.eval(knee->getMedialCondyle());
	}

    return result;
}


itk::Vector< double, 3 > ImplantsMatchFinalInfo::SetThicknessTibia(double medialThickness)
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

/*
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
*/

itk::Vector< double, 3 > ImplantsMatchFinalInfo::SetThicknessFemurCoronal(double medialThickness)
{
    double finalThickness = medialThickness;// -femurCartilage;
    Point movePoint = knee->getMedialCondyle();

    Plane planeTemp = femurImplant.getPosterior();
    planeTemp.reverseByPoint(femurImplant.getRodBasePoint());

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

/*
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
*/

itk::Vector< double, 3 > ImplantsMatchFinalInfo::SetThicknessFemurAxial(double medialThickness)
{
    double finalThickness = medialThickness;// -femurCartilage;
    Point movePoint = knee->getMedialInferiorFemurPoint();

    Plane planeTemp = femurImplant.getDistalPlane();
    planeTemp.reverseByPoint(femurImplant.getRodTopPoint());

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

/*
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
*/

ResectionThickness ImplantsMatchFinalInfo::GetFemurResectionAxial() const
{
    Plane planeTemp = femurImplant.getDistalPlane();
    planeTemp.reverseByPoint(femurImplant.getRodTopPoint(), false);
    Plane femur = TransformPlane(planeTemp, femurRotation, femurTranslation);

    ResectionThickness result;

	if (knee->getSurgerySide() == SurgerySideEnum::KLateral)
	{
		result.SurgerySide = femur.eval(knee->getLateralInferiorFemurPoint());
		result.NoSurgerySide = femur.eval(knee->getMedialInferiorFemurPoint());
	}
	else
	{
		result.NoSurgerySide = femur.eval(knee->getLateralInferiorFemurPoint());
		result.SurgerySide = femur.eval(knee->getMedialInferiorFemurPoint());
	}

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
    std::cout << "Tibia resection surgical: " << tibia.SurgerySide << " No surgical: " << tibia.NoSurgerySide << std::endl;
    std::cout << "Femur resection axial surgical: " << axialFemur.SurgerySide << " No surgical: " << axialFemur.NoSurgerySide << std::endl;
    std::cout << "Femur resection coronal surgical: " << coronalFemur.SurgerySide << " No surgical: " << coronalFemur.NoSurgerySide << std::endl;

    std::cout << "Translation femur: " << femurTranslation << std::endl;
    std::cout << "Translation tibia: " << tibiaTranslation << std::endl;

    Point onBone = GetFemurBoneTEALine().first - GetFemurBoneTEALine().second;
    Point onImplant = GetFemurImplantTEALine().first - GetFemurImplantTEALine().second;
    onBone.normalice();
    onImplant.normalice();

    std::cout << "femur bone tea vector: " << onBone << std::endl;
    std::cout << "femur implant tea vector: " << onImplant << std::endl;
}