#include "Knee.hpp"
#include "Types.hpp"
#include "implants_export.h"

class RegistrationPointsVTK;

struct IMPLANTS_EXPORT RegistrationPoints
{
    std::vector<PointTypeITK> points;
    RegistrationPoints(std::vector<PointTypeITK> pPoints);
    RegistrationPoints(std::vector<Point> pPoints);
};

class IMPLANTS_EXPORT FindRegistrationPoints
{
private:

    Knee mKnee;

    /*
    RegistrationPointsVTK * data;

    std::vector<RegistrationPoints> GetRegistrationPointsTibiaTemplate(bool isLeft = true);

    std::vector<cv::Point3d> GetRegistrationPointsTibiaTemplateLikeCV(bool isLeft = true);

    std::vector<cv::Point3d> GetRegistrationPointsFemurTemplateLikeCV(bool isLeft = true);

    double RegisterFemurTest(const std::vector<cv::Point3d>& templatePoints, const Point& vectorTea, const Point& vectorAp, const Point& center, double pScale = 1.0);

    double RegisterTibiaPoints();

    double RegisterTibiaLeft(const std::vector<cv::Point3d>& templatePoints, const Point& vectorTea, const Point& vectorAp, const Point& center, double sizeAp, double sizeTea);

    double RegisterTibiaRight(const std::vector<cv::Point3d>& templatePoints, const Point& vectorTea, const Point& vectorAp, const Point& center, double sizeAp, double sizeTea);
    */

public:
    FindRegistrationPoints(const Knee& pKnee);
    ~FindRegistrationPoints();
    std::vector<RegistrationPoints> GetRegistrationPointsFemur(std::vector<Point>& pCheckPoints, double& pError);
    std::vector<RegistrationPoints> GetRegistrationPointsTibia(std::vector<Point>& pCheckPoints, double& pError);
};