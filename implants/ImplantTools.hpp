#ifndef IMPLANT_TOOLS_H
#define IMPLANT_TOOLS_H

#include "Plane.hpp"
#include "Utils.hpp"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkImplicitPolyDataDistance.h"
#include <list>
#include <itkRigid3DTransform.h>

namespace ImplantTools
{
    struct Poly
    {
        std::vector<double> coeff;
        double minX, maxX;
        double Xc, Yc;
        double eval(double x);
        double Z;
        Point getLocalMin();
        Point getLocalMax();
        bool isFine = false;
    };

	struct PolyConstraintPoint
	{
		Point a, b;
		bool putConstraid = false;
	};

    cv::Mat GetRotateZ(const Point& vector);

    cv::Mat GetRotateY(const Point& vector);

    cv::Mat GetRotateX(const Point& vector);

    cv::Mat getRotateMatrix(const Point& axis, double angle);

    cv::Mat GetGeneralRotateTransformVectors(const Point pFromVector, const Point pToVector);

    double getAngleBetweenVectors(const Point& a, const Point& b);

    double getAngleBetweenVectorsDegree(const Point& a, const Point& b);

    double getDistanceBetweenPoints(const Point& a, const Point& b, bool square = false);

    double getDistanceBetweenPoints(const double a[3], const double b[3], bool square = false);

    //////////////////////////////////////////////////////////////////////it will remove

    Point getLocalMinimum(const std::vector<Point>& pPoints, const Plane& onPlane, const Point& vectorY, const Point& nearMin);

    Point getLocalMinimum(const std::vector<Point>& pPoints, const Plane& onPlane, const Point& vectorY, const Plane& nearMinPlane);

    Point getLocalMinimumTest(const std::list<Point>& pPoints, const Plane& onPlane, const Point& vectorY);
    //////////////////////////////////////////////////////////////////

    Point getLocalMinimum(const std::list<Point>& pPoints, const Plane& onPlane, const Point& vectorY, int degree = 7);

    Point getLocalMinimum(const std::vector<Point>& pPoints, const Plane& onPlane, const Point& vectorY, int degree = 7);

    Point getLocalMax(const std::list<Point>& pPoints, const Plane& onPlane, const Point& vectorY, int degree = 7);

    Point getLocalMax(const std::vector<Point>& pPoints, const Plane& onPlane, const Point& vectorY, int degree = 7);

    Poly polyFit(const std::list<Point>& pPoints, const cv::Mat& pTransformXY, int order, PolyConstraintPoint constraint);

    Poly polyFit(const std::vector<Point>& pPoints, const cv::Mat& pTransformXY, int order, PolyConstraintPoint constraint);

    vtkSmartPointer<vtkPolyData> getMaxContour(const vtkSmartPointer<vtkPolyData> polyData, const Point& pNormal, const Point& pPoint);

    vtkSmartPointer<vtkPolyData> getContours(const vtkSmartPointer<vtkPolyData> polyData, const Point& pNormal, const Point& pPoint);

	std::vector<std::pair<vtkSmartPointer<vtkPolyData>, vtkSmartPointer<vtkPoints>>> getAllContours(const vtkSmartPointer<vtkPolyData> polyData, const Point& pNormal, const Point& pPoint);

    void ExtractSortLines(const vtkSmartPointer<vtkPolyData> polyData, std::list<std::pair<vtkIdType, vtkIdType>>& lines);

    vtkIdType GetNearestPoints(const vtkSmartPointer<vtkPolyData> poly, const Point& pPoint);

    vtkIdType GetNearestPoints(const vtkSmartPointer<vtkPolyData> poly, const Line& pLine, const Plane& pPlane1, const Plane& pPlane2);

	std::pair<double, Point> GetDistancePlaneToSurface(const vtkSmartPointer<vtkPolyData> poly, const Plane& pPlane, const Plane& pUseOneSide = Plane());

    bool GetInterceptionWithSegment(const vtkSmartPointer<vtkPolyData> poly, const Point& p1, const Point& p2, Point& result);

    double GetInterceptionWithLine(const vtkSmartPointer<vtkImplicitPolyDataDistance> polyDistance, const Point& pRefPoint, const Point& p2, Point& result);

    Point GetFarestPoint(const vtkSmartPointer<vtkPolyData> poly, const Plane& pPlane);

    Point GetFarestPoint(const vtkSmartPointer<vtkPolyData> poly, const Plane& pPlane, const Plane& pCondition);

    void GetPointsOnContourSort(const vtkSmartPointer<vtkPolyData> contour, const std::vector<Plane>& pCondition, const Plane& nearPlane, std::vector<Point>& pResult);

    void GetPointsOnContourSort(const vtkSmartPointer<vtkPolyData> contour, const std::vector<Plane>& pCondition, const Point& nearPoint, std::vector<Point>& pResult);

    bool isPointInsideSphere(const Point& center, double radius, const Point& pPoint);

    void saveVectorPoints(const std::vector<Point>& pPoints, const std::string& pName, int axis = 3);

    std::pair<cv::Point2d, double> findCircle(const cv::Point2d& p1, const cv::Point2d& p2, const cv::Point2d& p3);

    std::pair<Point, double> fitSphere(const std::vector<Point>& pPoints);

    std::pair<double, double> fitEllipse(const vtkSmartPointer<vtkPolyData> pContour, const Point& pNormal, Point& center);

    int orientation(const cv::Point2d& p1, const cv::Point2d& p2, const cv::Point2d& p3); //0 collinear, 1 Clockwise, 2 Counterclockwise

    int orientation(const Point& p1, const Point& p2, const Point& p3, const cv::Mat& pRotationZ); //0 collinear, 1 Clockwise, 2 Counterclockwise
    
    int GetCornerPointOnContour(const std::vector<Point>& pPoints, const Point& pCenter, const Point& pExtremeA, const Point& pExtremeB, Plane pConditionPlane = Plane(), double pWeightA = 1, double pWeightB = 1);
    
    cv::Mat Rigid3DTransformToCV(const itk::Rigid3DTransform<>::Pointer transform);

    cv::Mat Rigid3DTransformToCVRotation(const itk::Rigid3DTransform<>::Pointer transform);

    cv::Mat Rigid3DTransformToCVTranslation(const itk::Rigid3DTransform<>::Pointer transform);

    cv::Mat JoinRigidTransform(const cv::Mat& rotation, const cv::Mat& translation);

    itk::Rigid3DTransform<>::Pointer getITKTransformFromCV(const cv::Mat& fullTransform);

    void squareCorner(const Point& pBasePoint, const Point& pCenter, const Point& pSecondPoint, std::vector<Point>& pData, double pSlopeAngle = 100);

    void squareCorner(const Point& pBasePoint, const Point& pCenter, const Point& pSecondPoint, const Plane& pRightSide, std::vector<Point>& pData, double pSlopeAngle = 100);

    Line GetSquareCornerFeatures(const Point& pBasePoint, const Point& pCenter, const Point& pSidePoint, const std::vector<Point>& pData, int& pDataSidePos, double pSlopeAngle = 100);

    void replaceVectorRange(const std::vector<Point>& pNewValues, int pInitPos, int pEndPos, std::vector<Point>& pVectorToReplace, bool sameOrder, bool includeInitPos = true);
    
    Point getHighestPointsOnTibia(const std::vector<Point>& pData, const Plane& pSplitLatPlane, Line& pTopLine, int& latPos, int& medPos);
    
    Plane TransformPlane(const Plane& plane, const itk::Rigid3DTransform<>::Pointer pTransform);

    std::vector<Point> getSortPointVectorFill(const std::vector<Point>& pVector, double pDistance);

	bool areBothSetOfPointsSeparated(const std::vector<Point>& pPoints1, const std::vector<Point>& pPoints2, const cv::Mat& myRotationZ, float margin = 0);
    
	std::vector<Point> increaseVectorPoints(const std::vector<Point>& pPoints, int beginPos, int endPos, float distance = 1.);

	void show(const vtkSmartPointer<vtkPolyData> poly1, const vtkSmartPointer<vtkPolyData> poly2);

    void show(vtkSmartPointer<vtkPolyData> poly, const std::vector<Point>& points, bool makePolyLine = false);

    vtkSmartPointer<vtkPolyData> getPolyLine(const std::vector<Point>& sortPoints);
};

#endif