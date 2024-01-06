#ifndef IMPLANTS_RIGISTRATION_POINTS
#define IMPLANTS_RIGISTRATION_POINTS

#include "Knee.hpp"
#include "Types.hpp"
#include "Utils.hpp"
#include "uka_implants_export.h"

namespace UKA
{
	namespace IMPLANTS
	{

		class RegistrationPointsVTK;

		struct UKA_IMPLANTS_EXPORT RegistrationPoints
		{
			std::vector<PointTypeITK> points;
			RegistrationPoints(std::vector<PointTypeITK> pPoints);
			RegistrationPoints(std::vector<Point> pPoints);
		};

		class UKA_IMPLANTS_EXPORT FindRegistrationPoints
		{
		private:
			Knee mKnee;

		public:
			FindRegistrationPoints(const Knee& pKnee);
			~FindRegistrationPoints();
			std::vector<RegistrationPoints> GetRegistrationPointsFemur(std::vector<Point>& pCheckPoints, double& pError);
			std::vector<RegistrationPoints> GetRegistrationPointsTibia(std::vector<Point>& pCheckPoints, double& pError);
		};
	}
}


#endif