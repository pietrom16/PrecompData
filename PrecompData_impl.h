/**  PrecompData_impl.h

	Copyright 2016 Pietro Mele
	Apache License 2.0
*/

#include "PrecompData.h"
#include <algorithm>
#include <cmath>
#include <functional>

namespace Utilities {

template<typename T, typename U, int nx, int ny>
PrecompData<T, U, nx, ny>::PrecompData()
	: interpolation(0), status(0), overSampling(2.0f)
{
}

template<typename T, typename U, int nx, int ny>
PrecompData<T, U, nx, ny>::PrecompData(const std::string _funcName)
	: interpolation(0), status(0), overSampling(2.0f), funcName(_funcName)
{
}


template<typename T, typename U, int nx, int ny>
int  PrecompData<T, U, nx, ny>::SetFunctionName(const std::string &_funcName)
{
	funcName = _funcName;
	return 0;
}

template<typename T, typename U, int nx, int ny>
int  PrecompData<T, U, nx, ny>::SetComment(const std::string &_comment)
{
	comment = _comment;
	return 0;
}

template<typename T, typename U, int nx, int ny>
std::string  PrecompData<T, U, nx, ny>::FunctionName() const
{
	return funcName;
}

template<typename T, typename U, int nx, int ny>
std::string  PrecompData<T, U, nx, ny>::Comment() const
{
	return comment;
}

template<typename T, typename U, int nx, int ny>
int  PrecompData<T, U, nx, ny>::SetOversampling(float ovs)
{
    if(ovs < 1.0f)
        return wrn_invalid_oversampling;

    overSampling = ovs;
    return 0;
}


// Precompute constant values

template<typename T, typename U, int nx, int ny>
int PrecompData<T, U, nx, ny>::PreComputeValues()
{
    //+CHECK
    // Set up conversion constants
    kRealInt = yData.size()/(xMax - xMin);
    kIntReal = 1/kRealInt;

    return 0;
}


// Coordinate <--> index transformations

template<typename T, typename U, int nx, int ny>
size_t PrecompData<T, U, nx, ny>::RtoI(T x) const     // real --> integer/index
{
    return size_t(kRealInt*(x - xMin));
}

template<typename T, typename U, int nx, int ny>
T PrecompData<T, U, nx, ny>::ItoR(size_t i) const     // integer/index --> real
{
    return xMin + kIntReal*T(i);
}


/// Data loading

// Regular grid, computed

template<typename T, typename U, int nx, int ny>
size_t  PrecompData<T, U, nx, ny>::Set(T (*Func1)(T x),
                                    T xmin, T xmax, size_t nPoints)     // line
{
	xMin = xmin;
	xMax = xmax;
	
    xData.resize(nPoints);
	yData.resize(nPoints);
    
	step = (xMax - xMin)/nPoints;

	T x = xMin;

	for(size_t i = 0; i < nPoints; ++i)
	{
		const T y = Func1(x);

        xData[i] = x;
		yData[i] = y;

		x += step;
	}
	
	PreComputeValues();

	return yData.size();
}


template<typename T, typename U, int nx, int ny>
size_t  PrecompData<T, U, nx, ny>::Set(T (*Func2)(T x, T y),
                                       T xmin, T xmax, size_t xnPoints,
                                       T ymin, T ymax, size_t ynPoints)    // plane
{} //+TODO

template<typename T, typename U, int nx, int ny>
size_t  PrecompData<T, U, nx, ny>::Set(T (*Func3)(T x, T y, T z),
                                    T xmin, T xmax, size_t xnPoints,
                                    T ymin, T ymax, size_t ynPoints,
                                    T zmin, T zmax, size_t znPoints)    // volume
{} //+TODO


/** AutoSet() : Automatic irregular grid, computed
 *
 *  Irregular grid:
 *    - Pick 25% of the points on a regular grid, for minimum uniform coverage.
 *    - Pick the remaining points from regions of the function with largest second derivative.
 */

template<typename T, typename U, int nx, int ny>
size_t  PrecompData<T, U, nx, ny>::AutoSet(T (*Func1)(T x), T xmin, T xmax, size_t nPoints)     // line
{
    xMin = xmin;
    xMax = xmax;

    PickBestPoints(Func1, nPoints, overSampling);

    PreComputeValues();

    return yData.size();
}


// Irregular grid, computed

template<typename T, typename U, int nx, int ny>
size_t  PrecompData<T, U, nx, ny>::Set(T (*Func1)(T x), const std::vector<T> &x)      // line
{} //+TODO

template<typename T, typename U, int nx, int ny>
size_t  PrecompData<T, U, nx, ny>::Set(T (*Func2)(T x, T y),
                                    const std::vector<T> &x,
                                    const std::vector<T> &y)      // plane
{} //+TODO

template<typename T, typename U, int nx, int ny>
size_t  PrecompData<T, U, nx, ny>::Set(T (*Func3)(T x, T y, T z),
                                    const std::vector<T> &x,
                                    const std::vector<T> &y,
                                    const std::vector<T> &z)      // volume
{} //+TODO

//+TODO Regular grid, load from file
//+TODO Irregular grid, load from file

/// Data retrieval

// Range UNchecked, 0 degree interpolation accessors

template<typename T, typename U, int nx, int ny>
T PrecompData<T, U, nx, ny>::operator()(T x) const
{
    return yData[RtoI(x)];
}

template<typename T, typename U, int nx, int ny>
T PrecompData<T, U, nx, ny>::operator()(T x, T y) const {} //+TODO

template<typename T, typename U, int nx, int ny>
T PrecompData<T, U, nx, ny>::operator()(T x, T y, T z) const {} //+TODO

// Range checked accessors; check Status()

template<typename T, typename U, int nx, int ny>
T PrecompData<T, U, nx, ny>::get(T x) {} //+TODO

template<typename T, typename U, int nx, int ny>
T PrecompData<T, U, nx, ny>::get(T x, T y) {} //+TODO

template<typename T, typename U, int nx, int ny>
T PrecompData<T, U, nx, ny>::get(T x, T y, T z) {} //+TODO


// Range checked accessors, interpolated; check Status()

template<typename T, typename U, int nx, int ny>
T PrecompData<T, U, nx, ny>::Interpolate(T x)
{
    RangeCheck(x);

    const size_t i = RtoI(x);

    if(i >= yData.size() - 2)
        return yData.back();

    const T x0 = ItoR(i);
    const T x1 = ItoR(i + 1);

    const T y = yData[i] + (yData[i + 1] - yData[i])*(x - x0)/(x1 - x0);

    return y;
}


template<typename T, typename U, int nx, int ny>
T PrecompData<T, U, nx, ny>::Interpolate(T x, T y) {} //+TODO

template<typename T, typename U, int nx, int ny>
T PrecompData<T, U, nx, ny>::Interpolate(T x, T y, T z) {} //+TODO


template<typename T, typename U, int nx, int ny>
void PrecompData<T, U, nx, ny>::Interpolation(int order)
{
    interpolation = order;
}


// Range check

template<typename T, typename U, int nx, int ny>
int PrecompData<T, U, nx, ny>::RangeCheck(T x)
{
    status = 0;

    if(x < xMin)
        status = wrn_x_less_than_min;
    else if(x > xMax)
        status = wrn_x_more_than_max;

    return status;
}


// Get the whole value set

template<typename T, typename U, int nx, int ny>
int PrecompData<T, U, nx, ny>::Get(std::vector<T> &_xData , std::vector<T> &_yData) const
{
    _xData = xData;
    _yData = yData;
    return 0;
}


/// Math functions

template<typename T, typename U, int nx, int ny>
T PrecompData<T, U, nx, ny>::FirstDerivative(T x1, T y1, T x2, T y2) const
{
    /// First derivative (central differences):  d1 = [f(x+1) - f(x)] / [(x+1) - x]
    return  (y2 - y1)/(x2 - x1);
}


template<typename T, typename U, int nx, int ny>
T PrecompData<T, U, nx, ny>::SecondDerivative(T x1, T y1, T x2, T y2, T x3, T y3) const
{
    /// Second derivative (central differences):  d2 = [f(x-1) - 2f(x) + f(x+1)] / {[(x+1) - (x-1)]/2}^2
    return  (y1 - 2*y2 + y3)/std::pow(0.5f*(x3 - x1), 2);
}


template<typename T, typename U, int nx, int ny>
int PrecompData<T, U, nx, ny>::PickBestPoints(T (*Func1)(T x), const size_t nPoints, const float overSampling)
{
    struct Point {
        T x, y;
        Point(T _x = 0.0, T _y = 0.0) : x(_x), y(_y) {}
        bool operator< (const Point &p) const { return x < p.x; }
        bool operator> (const Point &p) const { return x > p.x; }
    };

    struct PointCurv {     // abscissa and second derivative
        T x, d2;
        PointCurv(T _x = 0.0, T _d2 = 0.0) : x(_x), d2(_d2) {}
        bool operator> (const PointCurv &p) const { return fabs(d2) > fabs(p.d2); }
    };

    /// Oversample

    const size_t nSamples = size_t(overSampling*nPoints);
    const T step = (xMax - xMin)/nSamples;

    std::vector<PointCurv> samples;
    samples.reserve(nSamples);

    T x1 = xMin;
    T x2 = x1 + step;
    T x3 = x2 + step;
    T y1 = Func1(x1);
    T y2 = Func1(x2);
    T y3 = Func1(x3);

    PointCurv p;

    for(size_t i = 0; i < nSamples - 1; ++i)    // first and last points added later
    {
        p.x = x2;   // central point
        p.d2 = SecondDerivative(x1, y1, x2, y2, x3, y3);
        samples.push_back(p);
        
        x1 = x2;
        x2 = x3;
        x3 += step;
        y1 = y2;
        y2 = y3;
        y3 = Func1(x3);
    }

    // Sort based on decreasing second derivative absolute value
	std::sort(samples.begin(), samples.end(), std::greater<PointCurv>());

    /// Pick the points with highest second derivative (curvature)

    std::vector<Point> points;
    points.resize(nPoints);

    points[0].x = xMin;               // always get the first point

    for(size_t i = 1, j = 0; i < nPoints - 1; ++i, ++j)
        points[i].x = samples[j].x;

    points[nPoints - 1].x = xMax;     // always get the last point

    for(size_t i = 0; i < nPoints; ++i)
        points[i].y = Func1(points[i].x);

    // Sort the picked points based on increasing abscissa
    std::sort(points.begin(), points.end());

    // Copy data in member variables

    xData.clear();
    yData.clear();
    xData.resize(points.size());
    yData.resize(points.size());

    for(size_t i = 0; i < points.size(); ++i)
    {
        xData[i] = points[i].x;
        yData[i] = points[i].y;
    }

    return 0;
}


} // Utilities
