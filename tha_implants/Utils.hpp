#ifndef UTILS_H
#define UTILS_H

namespace THA
{
	namespace IMPLANTS
	{

		const double EPSILON = std::numeric_limits<double>::epsilon();
		const double PI = acos(-1.0);
		//const double LAT_MED_DIFFERENCE = 2.0;
		const double CEMENT = 1.0;

		struct FemurImplantInfo
		{
			double femurPosteriorLateralThickness;
			double femurPosteriorMedialThickness;
			double femurDistalLateralThickness;
			double femurDistalMedialThickness;
		};

		struct TibiaImplantInfo
		{
			double tibiaLateralThickness;
			double tibiaMedialThickness;
			//double slope;
		};
	}
}


#endif