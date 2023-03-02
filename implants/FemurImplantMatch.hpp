#ifndef FEMUR_IMPLANT_MATCH_H
#define FEMUR_IMPLANT_MATCH_H

#include <itkRigid3DTransform.h>
#include "FemurImplant.hpp"
#include "Knee.hpp"
#include "Plane.hpp"
#include "implants_export.h"
#include "vtkPolyData.h"

class IMPLANTS_EXPORT FemurImplantMatch
{
public:
	enum PlaneID{
		kPlaneA,//Femur posterior plane
		kPlaneB,//Femur posterior oblique plane
		kPlaneC,//Femur distal plane
		kPlaneD,//Femur anterior oblique plane
		kPlaneE,//Femur anterior plane
		kPlaneMid//Femur mid plane
	};

	FemurImplantMatch();

	void init(const FemurImplant& implant, const Knee& knee, bool useKneeCenterAlignment = true);

	Plane GetPlane(PlaneID id, bool translateByCondyle = true) const;

    itk::Matrix< double, 3, 3 > GetRotationMatrix() const;

    itk::Vector< double, 3 > GetTranslationMatrix() const;

    itk::Vector< double, 3 > GetTranslationMatrixByCortex() const;

	std::vector<Point> GetPointsNearPlane(PlaneID id, bool translateByCondyle = true, double distance = 0) const;

    std::vector<PointTypeITK> GetHullPoints(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut, PlaneID id, double distanceSide = 0, double distanceTop = 1.0, double angleLateral = 10, double angleMedial = 15, int amount = 200, bool posteriorLongCurve = false) const;

	Point TransformImplantPointToBone(const Point& pPoint, bool translateByCondyle = true) const;

    vtkSmartPointer<vtkPolyData> GetCuttingFemur(bool translateByCondyle);

protected:
	friend class Balance;

private:
	FemurImplant implant;
	Knee knee;
	cv::Mat rotationMatrix;
	cv::Mat translationMatrix;
	cv::Mat translationMatrixByCortex;
	bool isInit;
	bool useKneeCenterAlignment;

	void getRotationMatrix();

	bool getTranslationMatrix();

	bool getTranslationMatrixByCortex();

    vtkSmartPointer<vtkPolyData> getContour(const vtkSmartPointer<vtkPolyData> poly, const Point& pNormal, const Point& pPoint) const;

	Plane transformPlane(const Plane& plane, bool translateByCondyle = true) const;

    Plane finalTransformPlane(const Plane& plane, const itk::Rigid3DTransform<>::Pointer pTransform) const;

	Point getNearPointUnderCortex(const Plane& myPlane, std::vector<Point>& points, double distance = 0) const;

    void getVerticesCDE(const std::vector<Point>& points, const Point& downPoint, const Point& lateralPoint, const Point& medialPoint, const Point& topPoint, const Point& centerP1, const Point& centerP2, double distanceSide, double distanceTop, double angleLat, double angleMed, std::vector<Point>& vertices) const;

    void getVerticesA(const std::vector<Point>& points, const Point& downPoint, const Point& lateralPoint, const Point& medialPoint, const Point& topPoint, const Plane& midPlane, const Plane& currentPlane, const cv::Mat& pRotation, std::vector<Point>& vertices, double distance, int amount, bool longCurve) const;

    void getVerticesB(const std::vector<Point>& points, const Point& downPoint, const Point& lateralPoint, const Point& medialPoint, const Point& topPoint, const Plane& midPlane, const Plane& currentPlane, const cv::Mat& pRotation, std::vector<Point>& vertices, double distanceSide, double distanceTop, int amount) const;

    void separateMedialAndLateralPoints(const std::vector<Point>& points, const Plane& sagitalPlane, const Point& extremeLatPoint, const Point& extremeMedPoint, std::vector<Point>& lateralOut, std::vector<Point>& medialOut, Point& centerPointOut) const;
    
    bool getConcaveMainPoints(const std::vector<Point>& points, const Plane& currentPlane, const Plane& sagitalPlane, const Point& downPoint, int pDistanceTop, int pDistanceSideIn, std::vector<Point>& result, std::pair<int, int>& pInOutPos, double pMoveOnCenter = 0) const;

    std::pair<int, int> getSupExternalCornerPos(const std::vector<Point>& allPoints, const std::vector<Point>& contour, const Plane& currentPlane, const Plane& sagitalPlane, const Point& downPoint) const;
    
    std::vector<PointTypeITK> increaseVectorToAmount(const std::vector<Point>& points, int amount) const;

    Point movePointAtNormal(const Point& movePoint, const Point& nextPoint, const cv::Mat& rotationZ, double distance, bool changeMove = false, bool clockWise = false) const;

    Knee getKnee() const;
};


#endif
