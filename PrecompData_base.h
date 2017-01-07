/**  PrecompData_base.h

	Base class with utility members.
*/

#ifndef PRECOMPDATA_BASE_HPP
#define PRECOMPDATA_BASE_HPP

#include <string>

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
	int          SetOversampling(float ovs);
	float        Oversamping()  const { return overSampling; }

	int  Status() const { return status; }

	int  Interpolation() const { return interpolation; }
	void Interpolation(int order);

public:
	// Error return values
	static const size_t err_no_data              = -1,
	                    err_device_not_available = -2;

	// Warning return values
	static const size_t wrn_x_less_than_min      = -101,
	                    wrn_x_more_than_max      = -102,
	                    wrn_invalid_oversampling = -103;

public:
	// Test
	friend class PrecompData_test;

protected:

	std::string  funcName, comment;
	int          interpolation;
	int          status;

	bool    regularGrid;            // true if points are equally spaced on all axes
	float   overSampling;
};


} // Utilities


#endif // PRECOMPDATA_BASE_HPP
