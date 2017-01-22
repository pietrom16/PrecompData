/**  PrecompData.h

    Version 1.3.1

	Copyright 2016 Pietro Mele
	Apache License 2.0

	---
	
	Container of precalculated n-dimensional values, returned with interpolation/regression.
	Purpose: improve performance avoiding the realtime computation of complex functions.
	
	- Use constant tables for non linear functions' values.
	
    - Approximations:
        - Functions: Y = f(X),  where X and Y are vectors.
		- Interpolations: step, linear [TODO: cubic, spline, ...]
		- Regressions: [TODO: linear, minimal least squares, ...]
        - Extrapolations: [TODO: linear, polynomial, ...]
	
    - Store data in RAM and/or device memory.

    To enable device memory usage (e.g. GPGPU), define the PRECOMPDATA_DEVICE macro.
	
    Tradeoff: Templates do not allow to have objects with number of dimensions defined at run time.
              However they provide better run time performance.

	Warning: the function to be interpolated had better have a certain complexity,
	         otherwise it may be better to compute it straight away.
*/

#ifndef PRECOMP_DATA_H
#define PRECOMP_DATA_H

//#define PRECOMPDATA_DEVICE

#include "PrecompData_base.h"

#include <array>
#include <string>

#ifdef PRECOMPDATA_DEVICE
#include <vector>
#include <boost/compute/algorithm/copy.hpp>
#include <boost/compute/container/vector.hpp>
#endif

namespace Utilities {



/** PrecompData
	Set of points approximating a function with one independent variable
	an n dependent variables.
	f: x --> Y
  */

template<
    int nPoints = 100,     /* number of points to precompute  */
    typename TX = float,   /* data type of independent vector */
    typename TY = float,   /* data type of dependent vector   */
    int ny = 1             /* number of dimensions of the dependent vector */
>
class PrecompData : public PrecompData_base
{
public:

	// Data types for the dependent data and index
	typedef std::array<TY, ny> YData;

public:

	PrecompData();
	PrecompData(const std::string _funcName);

	int nxDimensions() const { return 1;  }
	int nyDimensions() const { return ny; }

    // Precompute constant values
	int PreComputeValues();

	// Coordinate <--> index transformation
	size_t ScalarToIndex(TX x)     const;     // scalar --> index
	TX     IndexToScalar(size_t i) const;     // index  --> scalar

	/// Data loading
	
	// Regular grid, computed
	size_t  set(YData (*Func)(TX x), TX xmin, TX xmax);
	size_t  set(TY (*Func)(TX x), TX xmin, TX xmax);

    // Automatic irregular grid, computed
	size_t  AutoSet(YData (*Func)(TX x), TX xmin, TX xmax);
	size_t  AutoSet(TY (*Func)(TX x), TX xmin, TX xmax);

    // Regular grid, load from file
	size_t  set(const std::string &dataFilename, TX xmin, TX xmax);

    // Irregular grid, load from file
	size_t  set(const std::string &dataFilename);      // grid contained in the file

	/// Data retrieval

	// Range UNchecked, 0 degree interpolation accessors
	size_t operator()(TX _x, TY &_y) const;
	size_t operator()(TX _x, YData &_y) const;
	TY     operator()(TX _x) const;

	// Range checked, 0 degree interpolation accessors
	size_t get(TX _x, TY &_y) const;
	size_t get(TX _x, YData &_y) const;

	// Range checked accessors, linear interpolation
	size_t Interpolate(TX _x, TY &_y) const;
	size_t Interpolate(TX _x, YData &_y) const;
	TY     Interpolate(TX _x);

	int RangeCheck(TX _x);

    // Get the whole value set
	int NPoints() const { return nPoints; }
	int get(std::vector<TX> &_xData, std::vector<TY> &_yData) const;
	int get(std::vector<TX> &_xData, std::vector<YData> &_yData) const;
	int Dump(int n = 0) const;
	int DumpElement(size_t j) const;

    // Evaluate error
	YData  EvaluateErrorKnownData() const;              // error on each dimension on known data
	TY     EvaluateAbsErrorKnownData() const;           // absolute error on known data
	YData  EvaluateError(int nTestPoints) const;        // error on each dimension on random points
	TY     EvaluateAbsError(int nTestPoints) const;     // absolute error on random points

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

protected:
	TY Norm(const YData&) const;
	TY FirstDerivative(TX x1, TY y1, TX x2, TY y2) const;
	TY SecondDerivative(TX x1, TY y1, TX x2, TY y2, TX x3, TY y3) const;
	TY SecondDerivativeAvg(TX x1, YData y1, TX x2, YData y2, TX x3, YData y3) const;
	TY SecondDerivativeAbsMax(TX x1, YData y1, TX x2, YData y2, TX x3, YData y3) const;
	int PickBestPoints(YData (*Func)(TX x), const float overSampling = 2.0f);
	int PickBestPoints(TY (*Func)(TX x), const float overSampling = 2.0f);

private:
	
	YData (*FuncX1Yn)(TX _x);
	TY    (*FuncX1Y1)(TX _x);

	std::array<TX, nPoints>     xData;
	std::array<YData, nPoints>  yData;

	TX  min, max, step;
	TX  kRealInt, kIntReal;     // conversion factors

#ifdef PRECOMPDATA_DEVICE
    boost::compute::vector<T>  device_line;
#endif

	friend class PrecompData_test;
};


} // Utilities


// Implementation include files
#include "PrecompData_impl.h"
#include "PrecompDataDevice_impl.h"
#include "PrecompData_dump.h"

#endif // PRECOMP_DATA_H
