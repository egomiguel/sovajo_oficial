#include "PatellaImplantMatchInfo.hpp"
#include "ImplantsException.hpp"
#include "ImplantTools.hpp"


PatellaImplantMatchInfo::PatellaImplantMatchInfo()
{
    isInit = false;
}

void PatellaImplantMatchInfo::init(const PatellaImplant& implant, const Knee& knee, const itk::Rigid3DTransform<>::Pointer pImplantToBonePatellaTransform)
{
    if (isInit == true)
    {
        throw ImplantExceptionCode::ALREADY_INITIALIZED_PATELLA_IMPLANT_MATCH_INFO;
    }
    this->implant = implant;
    this->knee = knee;

    patellaRotation = Rigid3DTransformToCVRotation(pImplantToBonePatellaTransform);
    patellaTranslation = Rigid3DTransformToCVTranslation(pImplantToBonePatellaTransform);

    isInit = true;
}

void PatellaImplantMatchInfo::setPatellaTransform(const itk::Rigid3DTransform<>::Pointer pImplantToBonePatellaTransform)
{
    patellaRotation = Rigid3DTransformToCVRotation(pImplantToBonePatellaTransform);
    patellaTranslation = Rigid3DTransformToCVTranslation(pImplantToBonePatellaTransform);
}

Plane PatellaImplantMatchInfo::TransformPlane(const Plane& plane, const cv::Mat& rotation, const cv::Mat& translation) const
{
    cv::Mat transformNormalVector = rotation * plane.getNormalVectorMat();
    cv::Mat transformPoint = (rotation * plane.getPointMat()) + translation;
    Plane transformPlane;
    transformPlane.init(Point(transformNormalVector), Point(transformPoint));
    return transformPlane;
}

cv::Mat PatellaImplantMatchInfo::Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer transform) const
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

cv::Mat PatellaImplantMatchInfo::Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer transform) const
{
    itk::Vector< double, 3 > translate = transform->GetOffset();
    cv::Mat result(3, 1, CV_64FC1);

    result.at <double>(0, 0) = translate[0];
    result.at <double>(1, 0) = translate[1];
    result.at <double>(2, 0) = translate[2];
    return result;
}

double PatellaImplantMatchInfo::getCutThickness() const
{
    Plane cutPlane = TransformPlane(implant.getBasePlane(), patellaRotation, patellaTranslation);
    return cutPlane.eval(knee.getPatellaDistalPosteriorPoint());
}

double PatellaImplantMatchInfo::getAngleML() const
{
    Plane cutPlane = TransformPlane(implant.getBasePlane(), patellaRotation, patellaTranslation);
    Plane frontPlane = knee.getPatellaFrontPlane();
    frontPlane.movePlane(knee.getPatellaCenter());

    Point vectorML = knee.getPatella().getPatellaLateralVector();

    Point vectorOficialLat = frontPlane.getProjectionVector(vectorML);
    Point vectorCutLat = cutPlane.getProjectionVector(vectorOficialLat);

    double angle = ImplantTools::getAngleBetweenVectorsDegree(vectorOficialLat, vectorCutLat);

    vectorCutLat.normalice();
    Point refPoint = knee.getPatellaCenter() + vectorCutLat;

    if (frontPlane.eval(refPoint) >= 0)
    {
        return angle;
    }
    else
    {
        return -angle;
    }
}

double PatellaImplantMatchInfo::getAngleSI() const
{
    Plane cutPlane = TransformPlane(implant.getBasePlane(), patellaRotation, patellaTranslation);
    Plane frontPlane = knee.getPatellaFrontPlane();
    frontPlane.movePlane(knee.getPatellaCenter());

    Point vectorSI = knee.getPatella().getPatellaInferiorVector();

    Point vectorOficialInf = frontPlane.getProjectionVector(vectorSI);
    Point vectorCutInf = cutPlane.getProjectionVector(vectorOficialInf);

    double angle = ImplantTools::getAngleBetweenVectorsDegree(vectorOficialInf, vectorCutInf);
    
    vectorCutInf.normalice();
    Point refPoint = knee.getPatellaCenter() + vectorCutInf;

    if (frontPlane.eval(refPoint) >= 0)
    {
        return angle;
    }
    else
    {
        return -angle;
    }
}

