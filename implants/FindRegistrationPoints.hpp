#ifndef IMPLANTS_RIGISTRATION_POINTS
#define IMPLANTS_RIGISTRATION_POINTS

#include "Knee.hpp"
#include "Types.hpp"
#include "tka_implants_export.h"

namespace TKA
{
	namespace IMPLANTS
	{

		class RegistrationPointsVTK;

		struct TKA_IMPLANTS_EXPORT RegistrationPoints
		{
			std::vector<PointTypeITK> points;
			RegistrationPoints(std::vector<PointTypeITK> pPoints);
			RegistrationPoints(std::vector<Point> pPoints);
		};

		class TKA_IMPLANTS_EXPORT FindRegistrationPoints
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
	}
}


#endif