/**  PrecompData.h

    Version 1.2.0

	Copyright 2016 Pietro Mele
	Apache License 2.0

	---
	
	Container of precalculated values, returned with interpolation/regression.
	Purpose: improve performance avoiding the realtime computation of complex functions.
	
	- Use constant tables for non linear functions' values.
	- Approximations:
		- Interpolations: linear, cubic, spline, ...
		- Regressions: linear, cubic, ...
	- Store data in RAM and/or device memory.

    To enable device memory usage, define the PRECOMPDATA_DEVICE macro.
	
	Warning: the function to be interpolated had better have a certain complexity,
	         otherwise it may be better to compute it straight away.
*/

#ifndef PRECOMP_DATA_H
#define PRECOMP_DATA_H

//#define PRECOMPDATA_DEVICE

#include <string>
#include <vector>

#ifdef PRECOMPDATA_DEVICE
#include <boost/compute/algorithm/copy.hpp>
#include <boost/compute/container/vector.hpp>
#endif

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

    // Automatic irregular grid, computed
    size_t  AutoSet(T (*Func1)(T x), T xmin, T xmax, size_t nPoints = 100);     // line: 1D --> 1D

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
	void Interpolation(int order);

    int RangeCheck(T x);

    // Get the whole value set
    int Get(std::vector<T> &_xData, std::vector<T> &_yData) const;

    /// GPGPU

#ifdef PRECOMPDATA_DEVICE
    int InitDevice();

    int CopyOnDevice(boost::compute::device         &device,
                     boost::compute::context        &context,
                     boost::compute::command_queue  &queue,
                     boost::compute::vector<T>      &device_line);

    int CopyOnDevice(boost::compute::device         **device      = nullptr,
                     boost::compute::context        **context     = nullptr,
                     boost::compute::command_queue  **queue       = nullptr,
                     boost::compute::vector<T>      **device_line = nullptr);

    // Copy a subset
	int CopyOnDevice(T xbeg, T xend);
	int CopyOnDevice(T xbeg, T xend, T ybeg, T yend);
	int CopyOnDevice(T xbeg, T xend, T ybeg, T yend, T zbeg, T zend);
#endif // PRECOMPDATA_DEVICE

public:
    // Error return values
    static const int err_no_data              = -1,
                     err_device_not_available = -2;
    
    // Warning return values
    static const int wrn_x_less_than_min      = -101,
                     wrn_x_more_than_max      = -102;

protected:
    T FirstDerivative(T x1, T y1, T x2, T y2) const;
    T SecondDerivative(T x1, T y1, T x2, T y2, T x3, T y3) const;
    int PickBestPoints(T (*Func1)(T x), const size_t nPoints, const int overSampling = 2);

private:
	
	std::string     funcName, comment;
	int             interpolation;
	int             status;

    std::vector<T>  xData;
    std::vector<T>  yData;
    T               xMin, xMax, step;
	T               kRealInt, kIntReal;     // conversion factors


	//+TODO: plane, volume

#ifdef PRECOMPDATA_DEVICE
    boost::compute::vector<T>  device_line;
#endif

};

} // Utilities


// Implementation include files
#include "PrecompData_impl.h"
#include "PrecompDataDevice_impl.h"

#endif // PRECOMP_DATA_H
