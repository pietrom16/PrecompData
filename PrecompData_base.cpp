/**  PrecompData_base.cpp

	Base class with utility members.
*/

#include "PrecompData_base.h"


namespace Utilities {


const std::string PrecompData_base::version = PrecompData_VERSION_MAJOR "." PrecompData_VERSION_MINOR "." PrecompData_VERSION_PATCH;


PrecompData_base::PrecompData_base()
    : interpolation(0), status(0), regularGrid(false), overSampling(2.0f)
{
}


PrecompData_base::PrecompData_base(const std::string _funcName)
    : funcName(_funcName), interpolation(0), status(0), regularGrid(false), overSampling(2.0f)
{
}


int  PrecompData_base::SetFunctionName(const std::string &_funcName)
{
	funcName = _funcName;
	return 0;
}


int  PrecompData_base::SetComment(const std::string &_comment)
{
	comment = _comment;
	return 0;
}


std::string  PrecompData_base::FunctionName() const
{
	return funcName;
}


std::string  PrecompData_base::Comment() const
{
	return comment;
}


size_t  PrecompData_base::SetOversampling(float ovs)
{
	if(ovs < 1.0f)
		return wrn_invalid_oversampling;

	overSampling = ovs;
	return 0;
}


void PrecompData_base::Interpolation(int order)
{
	interpolation = order;
}



} // Utilities
