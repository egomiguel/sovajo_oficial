#ifndef UTILS_PKA_H
#define UTILS_PKA_H

namespace PKA
{
	namespace IMPLANTS
	{
		const double EPSILON = std::numeric_limits<double>::epsilon();
		const double PI = acos(-1.0);
		//const double LAT_MED_DIFFERENCE = 2.0;
		const double CEMENT = 1.0;

		struct FemurImplantInfo
		{
			double femurPosteriorThickness;
			double femurDistalThickness;
		};

		struct TibiaImplantInfo
		{
			double tibiaThickness;
			//double slope;
		};
	}
}

#endif