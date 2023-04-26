#ifndef PATELLA_IMPLANT_MATCH_H
#define PATELLA_IMPLANT_MATCH_H

#include "PatellaImplant.hpp"
#include "Knee.hpp"
#include "Plane.hpp"
#include "implants_export.h"

class IMPLANTS_EXPORT PatellaImplantMatch
{
public:
    PatellaImplantMatch();

	void init(const PatellaImplant& implant, const Knee& knee);

	Plane getPatellaCutPlane() const;

	Point transformImplantPoint(const Point& point) const;

    itk::Matrix< double, 3, 3 > GetRotationMatrix() const;

    itk::Vector< double, 3 > GetTranslationMatrix() const;

    std::vector<PointTypeITK> GetHullPoints(const itk::Rigid3DTransform<>::Pointer pTransformIn, itk::Rigid3DTransform<>::Pointer pTransformOut, double distance = 1.0, int amount = 200) const;

    vtkSmartPointer<vtkPolyData> GetCuttingPatella() const;

private:
    PatellaImplant implant;
	Knee knee;
	cv::Mat rotationMatrix;
	cv::Mat translationMatrix;
	bool isInit;

	void makeRotationMatrix();

	void makeTranslationMatrix();

	Plane transformPlane(const Plane& plane) const;

    Plane transformPlane(const Plane& plane, const itk::Rigid3DTransform<>::Pointer pTransform) const;

    Point movePointAtNormal(const Point& movePoint, const Point& nextPoint, const cv::Mat& rotationZ, double distance, bool changeMove = false, bool clockWise = true) const;

    std::vector<PointTypeITK> increaseVectorToAmount(const std::vector<Point>& points, int amount) const;
};


#endif
