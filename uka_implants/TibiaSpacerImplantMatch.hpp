#ifndef TIBIA_SPACER_IMPLANT_UKA_MATCH_H
#define TIBIA_SPACER_IMPLANT_UKA_MATCH_H

#include "TibiaImplant.hpp"
#include "TibiaSpacerImplant.hpp"
#include "Plane.hpp"
#include "uka_implants_export.h"

namespace UKA
{
	namespace IMPLANTS
	{

		class UKA_IMPLANTS_EXPORT TibiaSpacerImplantMatch
		{
		public:
			TibiaSpacerImplantMatch();

			void init(const TibiaImplant& tibiaImplant, const TibiaSpacerImplant& spacerImplant);

			itk::Matrix< double, 3, 3 > GetRotationMatrix() const;

			itk::Vector< double, 3 > GetTranslationMatrix() const;

		private:
			TibiaImplant implant;
			TibiaSpacerImplant spacerImplant;
			cv::Mat rotationMatrix;
			cv::Mat translationMatrix;
			bool isInit;

			void makeRotationMatrix();

			void makeTranslationMatrix();
		};
	}
}

#endif
