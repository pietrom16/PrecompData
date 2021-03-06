/**  PrecompData_multi.h

	Copyright 2016 Pietro Mele
	Apache License 2.0
*/

#ifndef PRECOMP_DATA_MULTI_H
#define PRECOMP_DATA_MULTI_H

//#define PRECOMPDATA_DEVICE

#include <string>
#include <boost/multi_array.hpp>

#ifdef PRECOMPDATA_DEVICE
#include <vector>
#include <boost/compute/algorithm/copy.hpp>
#include <boost/compute/container/vector.hpp>
#endif

namespace Utilities {



/** PrecompData
	Set of points approximating a multidimensional function/hypersurface
	f: X --> Y
  */

template<
    typename TX = float,   /* data type of independent vector */
    typename TY = float,   /* data type of dependent vector   */
    int nx,                /* number of dimensions of the independent vector */
    int ny                 /* number of dimensions of the dependent vector   */
>
class PrecompData
{
public:

	// Data types for the dependent and independent data and indices
	typedef boost::multi_array<TX, nx> XData;
	typedef boost::multi_array<TY, ny> YData;
	typedef typename XData::index X;
	typedef typename YData::index Y;

public:

	PrecompData();
	PrecompData(const std::string _funcName);

	int          SetFunctionName(const std::string &_funcName);
	int          SetComment(const std::string &_comment);
	std::string  FunctionName() const;
	std::string  Comment()      const;
	int          SetOversampling(float ovs);
	float        Oversamping()  const { return overSampling; }

	int          nxDimensions() const { return nx; }
	int          nyDimensions() const { return ny; }

	// Precompute constant values
	int PreComputeValues();

	// Coordinate <--> index transformation
	size_t VectorToIndex(X x)      const;     // vector --> index
	size_t ScalarToIndex(TX x)     const;     // scalar --> index
	X      IndexToVector(size_t i) const;     // index  --> vector
	TX     IndexToScalar(size_t i) const;     // index  --> scalar

	/// Data loading

	// Regular grid, computed
	size_t  Set(Y (*Func)(X x), X xmin, X xmax, size_t nPoints);
	size_t  Set(TY (*Func)(TX x), TX xmin, TX xmax, size_t nPoints);

	// Automatic irregular grid, computed
	size_t  AutoSet(Y (*Func)(X x), X xmin, X xmax, size_t nPoints = 100);
	size_t  AutoSet(TY (*Func)(TX x), TX xmin, TX xmax, size_t nPoints = 100);

	// Regular grid, load from file
	size_t  Set(const std::string &dataFilename, X xmin, X xmax);

	// Irregular grid, load from file
	size_t  Set(const std::string &dataFilename);      // grid contained in the file


	/// Data retrieval

	// Range UNchecked, 0 degree interpolation accessors
	Y  operator()(X x)  const;
	TY operator()(TX x) const;

	// Range checked accessors; check Status()
	Y get(X x);

	// Range checked accessors, interpolated; check Status()
	Y  Interpolate(X x);
	TY Interpolate(TX x);

	int Status() const { return status; }

	int  Interpolation() const { return interpolation; }
	void Interpolation(int order);

	int RangeCheck(X x);
	int RangeCheck(TX x);

	// Get the whole value set
	int Get(std::vector<X> &_xData, std::vector<Y> &_yData) const;
	int Get(std::vector<TX> &_xData, std::vector<TY> &_yData) const;
	int Dump(int n = 0) const;
	int DumpElement(size_t j) const;

	// Evaluate error
	Y  EvaluateErrorKnownData() const;              // error on each dimension on known data
	TY EvaluateAbsErrorKnownData() const;           // absolute error on known data
	Y  EvaluateError(int nTestPoints) const;        // error on each dimension on random points
	TY EvaluateAbsError(int nTestPoints) const;     // absolute error on random points

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
	                 wrn_x_more_than_max      = -102,
	                 wrn_invalid_oversampling = -103;

public:
	// Test
	friend class PrecompData_test;

protected:
	TY Norm(const Y&) const;
	TY FirstDerivative(TX x1, TY y1, TX x2, TY y2) const;
	TY SecondDerivative(TX x1, TY y1, TX x2, TY y2, TX x3, TY y3) const;
	int PickBestPoints(Y (*Func)(X x), const size_t nPoints, const float overSampling = 2.0f);
	int PickBestPoints(TY (*Func)(TX x), const size_t nPoints, const float overSampling = 2.0f);

private:

	std::string     funcName, comment;
	int             interpolation;
	int             status;

	Y  (*FuncX)(X x);
	TY (*FuncTX)(TX x);

	XData   xData;
	YData   yData;
	X       min, max, step;
	X       kRealInt, kIntReal;     // conversion factors

	bool    regularGrid;            // true if points are equally spaced on all axes
	float   overSampling;

#ifdef PRECOMPDATA_DEVICE
	boost::compute::vector<T>  device_line;
#endif

};


} // Utilities


// Implementation include files
#include "PrecompData_impl.h"
#include "PrecompDataDevice_impl.h"
#include "PrecompData_dump.h"

#endif // PRECOMP_DATA_MULTI_H
