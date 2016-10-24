/**  PrecompData.h

	Copyright 2016 Pietro Mele
	Apache License 2.0

	---
	
	Container of precalculated values, returned with interpolation/regression.
	Purpose: improve performance avoiding the realtime computation of complex functions.
	
	- Use constant tables for non linear functions' values.
	- Approximations:
		- Interpolations: linear, cubic, spline, ...
		- Regressions: linear, cubic, ...
	- Store data in RAM and/or GPU memory.
	
	Warning: the function to be interpolated had better have a certain complexity,
	         otherwise it may be better to compute it straight away.
*/

#ifndef PRECOMP_DATA_H
#define PRECOMP_DATA_H

#include <string>
#include <vector>

namespace Utilities {

template<typename T>
class PrecompData
{
public:

	PrecompData();
	PrecompData(const std::string _funcName);

	int          SetFunctionName(const std::string &_funcName);
	int          SetComment(const std::string &_comment);
	std::string  FunctionName() const;
	std::string  Comment()      const;
	
    // Precompute constant values
	int PreComputeValues();

	// Coordinate <--> index transformation
	size_t RtoI(T x)      const;     // real --> integer
	T      ItoR(size_t i) const;     // integer --> real
	
	/// Data loading
	
	// Regular grid, computed
    size_t  Set(T (*Func1)(T x),           T xmin, T xmax, size_t nPoints);     // line: 1D --> 1D
    size_t  Set(T (*Func2)(T x, T y),      T xmin, T xmax, size_t xnPoints,
	                                       T ymin, T ymax, size_t ynPoints);    // plane: 2D --> 1D
    size_t  Set(T (*Func3)(T x, T y, T z), T xmin, T xmax, size_t xnPoints,
	                                       T ymin, T ymax, size_t ynPoints,
	                                       T zmin, T zmax, size_t znPoints);    // volume: 3D --> 1D
	
	// Regular grid, load from file
    size_t  Set(const std::string &dataFilename, T xmin, T xmax);    // line
    size_t  Set(const std::string &dataFilename, T xmin, T xmax,
	                                             T ymin, T ymax);    // line
    size_t  Set(const std::string &dataFilename, T xmin, T xmax,
	                                             T ymin, T ymax,
	                                             T zmin, T zmax);    // volume

	// Irregular grid, computed
    size_t  Set(T (*Func1)(T x),           const std::vector<T> &x);    // line
    size_t  Set(T (*Func2)(T x, T y),      const std::vector<T> &x,
	                                       const std::vector<T> &y);    // plane
    size_t  Set(T (*Func3)(T x, T y, T z), const std::vector<T> &x,
	                                       const std::vector<T> &y,
	                                       const std::vector<T> &z);    // volume

	// Irregular grid, load from file
    size_t  Set(const std::string &dataFilename);      // grid contained in the file


	/// Data retrieval

	// Range UNchecked, 0 degree interpolation accessors
	T operator()(T x)           const;
	T operator()(T x, T y)      const;
	T operator()(T x, T y, T z) const;

	// Range checked accessors; check Status()
	T get(T x);
	T get(T x, T y);
	T get(T x, T y, T z);

	// Range checked accessors, interpolated; check Status()
	T Interpolate(T x);
	T Interpolate(T x, T y);
	T Interpolate(T x, T y, T z);

	int Status() const { return status; }

	int  Interpolation() const { return interpolation; }
	void Interpolate(int order = 0);

	/// GPGPU

	int CopyOnGPU();
	
	// Copy a subset
	int CopyOnGPU(T xbeg, T xend);
	int CopyOnGPU(T xbeg, T xend, T ybeg, T yend);
	int CopyOnGPU(T xbeg, T xend, T ybeg, T yend, T zbeg, T zend);
	
private:
	
	std::string     funcName, comment;
	int             interpolation;
	int             status;

	std::vector<T>  line;
	T               xMin, xMax, step;
	T               kRealInt, kIntReal;     // conversion factors
	//+TODO: plane, volume
};

} // Utilities


#include "PrecompData_impl.h"	// implementation

#endif // PRECOMP_DATA_H
