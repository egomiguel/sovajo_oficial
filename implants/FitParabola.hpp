#include "Point.hpp"
#include "Plane.hpp"

class FitParabola
{
private:
    double a, b, c;
    std::vector<Point> points;
    cv::Mat zRotation;
    bool makeFit();
    double zCoord;
public:
    FitParabola(const std::vector<Point>& pPoints, const Plane& plane, bool& result);
    Point getVertex() const;
    std::vector<Point> getPoints(double rangeInit = -100, double rangeEnd = 100) const;
    static std::vector<Point> SmoothCurve(const std::vector<Point>& pPoints, int kernelSize);
};