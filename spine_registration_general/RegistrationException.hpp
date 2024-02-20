#ifndef EXCEPTION_REGISTRATION_H
#define EXCEPTION_REGISTRATION_H

#include <exception>
#include <string>
#include "spine_registration_export.h"

namespace SPINE
{
	namespace REGISTRATION
	{

		class SPINE_REGISTRATION_EXPORT RegistrationException : public std::exception
		{
			std::string msg_;
		public:
			RegistrationException(const std::string& msg) : msg_(msg) {}

			virtual const char* what() const noexcept override
			{
				return msg_.c_str();
			}
		};
	}
}

#endif