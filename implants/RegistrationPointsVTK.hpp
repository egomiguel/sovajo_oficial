
#include "Knee.hpp"

class RegistrationPointsVTK
{
private:
    Point lastMedialBorderFemur;
    Point mMedialTibiaRef, mLateralTibiaRef;
    bool mGoodMediaRef, mGoodLateralRef;
    Point myTuber;
    Point mTransverseRef, mObliqueRef, mTuberUpRef, mTuberDownRef, mTuberSideRef;

    Point mVectorFemurAxisUp, mVectorFemurAPFront, mVectorFemurTEALat;
    Point mVectorTibiaAxisUp, mVectorTibiaAPFront, mVectorTibiaTEALat;

    Point mBigCenterLat, mBigCenterMed, mCenterEllipse;
    Plane mTibiaTransverse;

    /*Point medialTibiaPoint1, lateralTibiaPoint1;
    Point medialTibiaPoint2, lateralTibiaPoint2;
    bool goodTibiaPoints;*/
    struct Sphere
    {
        Point center;
        double radius;
        bool real;
        Sphere(const Point& pCenter, double pRadius);
        Sphere();
        bool isPointInside(const Point& pnt) const;
    };
    Knee mKnee;
    vtkSmartPointer<vtkPolyData> medialSideFemur, lateralSideFemur, tibiaSidePoly;
    std::vector<Point> reduceSortPointsByRadius(const std::vector<Point>& sortPoints, double radius);
    void reduceSortPointsByRadius(const std::vector<Point>& sortPoints, double radius, std::vector<Point>& result);
    std::vector<Point> reduceSortPointsByAmount(const std::vector<Point>& sortPoints, int pointsAmount);
    Point getNearPointToLine(const Line& pLine, const vtkSmartPointer<vtkPoints> Ppoints, const Sphere& pExclude);
    std::vector<Point> GetMedialMainPointsFemur(int amount = 6, int pTest = 0);
    std::vector<Point> GetLateralMainPointsFemur(int amount = 6, int pTest = 0);
    //std::vector<Point> GetUpPointsTibia(const std::vector<Point>& pPoints, const Point& pPlateau);
    std::vector<Point> GetObliquePointsTibia(const Point& fixPoint, const Point& movePoint, int amount = 4);
    std::vector<Point> GetTransversalPointsTibia(const Point& fixPoint, const Point& movePoint, int amount = 4);
    void makeAxis();
    
public:
    RegistrationPointsVTK(const Knee& pKnee);
    std::vector<Point> GetLateralBorderPointsFemur(int amount = 6);
    std::vector<Point> GetMedialBorderPointsFemur(int amount = 6);
    std::vector<Point> GetLateralObliquePointsFemur(int amount = 3);
    std::vector<Point> GetSagitalPointsFemur(int amount = 3);

    /*std::pair<std::vector<Point>, std::vector<Point>> GetTransversePointsTibia(int amount = 4);
    std::pair<std::vector<Point>, std::vector<Point>> GetUpsPointsTibia(int amount = 4);
    std::vector<Point> GetMedialObliqueOut(int amount = 4);
    std::vector<Point> GetMedialObliqueIn(int amount = 3);
    std::vector<Point> GetLateralObliqueOut(int amount = 4);
    std::vector<Point> GetLateralObliqueIn(int amount = 3);*/

    std::pair<Point, double> getMedialPlateauCirclePoints(std::vector<Point>& circle, std::vector<Point>& transversePoints);
    std::pair<Point, double> getLateralPlateauCirclePoints(std::vector<Point>& circle);
    std::vector<Point> GetTibiaMedialObliquePoints(int amount = 4);
    std::vector<Point> GetTibiaLateralObliquePoints(int amount = 4);
    std::vector<Point> GetTibiaOnTuberAxisPoints(int amount = 4);
    std::vector<Point> GetTibiaOnTuberSidePoints(int amount = 4);

    /////////////////////////////////////////////////////////////////////

    std::vector<Point> GetTibiaTemplatePoints();
    std::vector<Point> GetTibiaTemplatePointsRight();
    Point GetTibiaCenterEllipse();

    void GetFemurObliqueLateralEpi(const std::vector<Point>& pRef, std::vector<Point>& result, double pRadius = 15);
    void GetFemurObliqueMedialEpi(const std::vector<Point>& pRef, std::vector<Point>& result, double pRadius = 15);
    std::vector<Point> GetFemurTemplatePoints(bool isLeft);
};