#ifndef EXCEPTION_ROD_H
#define EXCEPTION_ROD_H

#include <exception>
#include <string>
#include "CheckImageRod_export.h"

class CHECKIMAGEROD_EXPORT CheckImageRodException : public std::exception
{
    std::string msg_;
public:
    CheckImageRodException(const std::string& msg) : msg_(msg) {}

    virtual const char* what() const noexcept override
    {
        return msg_.c_str();
    }
};

#endif
