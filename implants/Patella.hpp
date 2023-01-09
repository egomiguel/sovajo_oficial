#ifndef PATELLA_H
#define PATELLA_H

#include "Plane.hpp"
#include "Utils.hpp"
#include <vector>
#include "implants_export.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkImplicitPolyDataDistance.h"

class IMPLANTS_EXPORT Patella
{
public:
    Patella();

    void init(const Point& pRotationPoint, const Point& pLateralPoint, const Point& pMedialPoint, const Point& pInferiorPoint, const vtkSmartPointer<vtkPolyData> pPatellaPoly );

    Plane getPatellaPlane(bool isRight, Point& furthestPoint) const;

    Point getPatellaCenter() const;

    Point getPatellaRotationPoint() const;

    Point getPatellaInferiorVector() const;

    Point getPatellaLateralVector() const;

    double getPatellaThickness() const;

    double getPatellaDiameter() const;

    bool getIsInit();

    vtkSmartPointer<vtkPolyData> getPatellaPoly() const;

private:
    //std::vector<Point> mPatellaPoints;
    vtkSmartPointer<vtkPolyData> mPatellaPoly;
    Point mRotationPoint, mLateralPoint, mMedialPoint, mInferiorPoint;
    bool isInit;
    Plane getFrontPlane(bool isRight) const;
    void getPatellaFront(const Point& a, const Point& b, const Point& center, const vtkSmartPointer<vtkImplicitPolyDataDistance>& polyDistance, std::vector<Point>& points) const;
};


#endif
