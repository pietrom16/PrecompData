/**  PrecompData_base.h

	Base class with utility members.
*/

#ifndef PRECOMPDATA_BASE_HPP
#define PRECOMPDATA_BASE_HPP

#include <string>

// Version
#define PrecompData_VERSION_MAJOR "1"
#define PrecompData_VERSION_MINOR "4"
#define PrecompData_VERSION_PATCH "0"
#define PrecompData_VERSION "1.4.0"

namespace Utilities {


class PrecompData_base
{
public:

	PrecompData_base();
	PrecompData_base(const std::string _funcName);

	int          SetFunctionName(const std::string &_funcName);
	int          SetComment(const std::string &_comment);
	std::string  FunctionName() const;
	std::string  Comment()      const;
	size_t       SetOversampling(float ovs);
	float        Oversamping()  const { return overSampling; }

	int  Interpolation() const { return interpolation; }
	void Interpolation(int order);

public:
	int  Status() const { return status; }

	// Error return values
	static const int err_no_data              = -1,
	                 err_no_function          = -2,
	                 err_device_not_available = -3;

	// Warning return values
	static const int wrn_x_less_than_min      = -101,
	                 wrn_x_more_than_max      = -102,
	                 wrn_invalid_oversampling = -103;

public:
	// Test
	friend class PrecompData_test;

	static const std::string version;

protected:

	std::string  funcName, comment;
	int          interpolation;
	int          status;

	bool    regularGrid;            // true if points are equally spaced on all axes
	float   overSampling;
};


} // Utilities


#endif // PRECOMPDATA_BASE_HPP
