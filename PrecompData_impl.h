/**  PrecompData_impl.h

	Copyright 2016 Pietro Mele
	Apache License 2.0
*/

#include <algorithm>
#include <cmath>
#include <functional>

namespace Utilities {

template<typename TX, typename TY, int nx, int ny>
PrecompData<TX, TY, nx, ny>::PrecompData()
	: interpolation(0), status(0), overSampling(2.0f)
{
}

template<typename TX, typename TY, int nx, int ny>
PrecompData<TX, TY, nx, ny>::PrecompData(const std::string _funcName)
	: interpolation(0), status(0), overSampling(2.0f), funcName(_funcName)
{
}


template<typename TX, typename TY, int nx, int ny>
int  PrecompData<TX, TY, nx, ny>::SetFunctionName(const std::string &_funcName)
{
	funcName = _funcName;
	return 0;
}

template<typename TX, typename TY, int nx, int ny>
int  PrecompData<TX, TY, nx, ny>::SetComment(const std::string &_comment)
{
	comment = _comment;
	return 0;
}

template<typename TX, typename TY, int nx, int ny>
std::string  PrecompData<TX, TY, nx, ny>::FunctionName() const
{
	return funcName;
}

template<typename TX, typename TY, int nx, int ny>
std::string  PrecompData<TX, TY, nx, ny>::Comment() const
{
	return comment;
}

template<typename TX, typename TY, int nx, int ny>
int  PrecompData<TX, TY, nx, ny>::SetOversampling(float ovs)
{
    if(ovs < 1.0f)
        return wrn_invalid_oversampling;

    overSampling = ovs;
    return 0;
}


// Precompute constant values

template<typename TX, typename TY, int nx, int ny>
int PrecompData<TX, TY, nx, ny>::PreComputeValues()
{
    //+CHECK
    // Set up conversion constants
    kRealInt = yData.size()/(xMax - xMin);
    kIntReal = 1/kRealInt;

    return 0;
}


// Coordinate <--> index transformations

template<typename TX, typename TY, int nx, int ny>
size_t PrecompData<TX, TY, nx, ny>::RtoI(X x) const     // real --> integer/index
{
    return size_t(kRealInt*(x - xMin));
}

/* //+?
template<typename TX, typename TY, int nx, int ny>
PrecompData<TX, TY, nx, ny>::X PrecompData<TX, TY, nx, ny>::ItoR(size_t i) const     // integer/index --> real
{
    return xMin + kIntReal*T(i);
}
*/
/// Data loading

// Regular grid, computed

template<typename TX, typename TY, int nx, int ny>
size_t  PrecompData<TX, TY, nx, ny>::Set(Y       (*Func)(X x), 
                                         X       xmin,
                                         X       xmax,
                                         size_t  nPoints)
{
    //+TODO
    min = xmin;
    max = xmax;

    xData.resize(nPoints);
    yData.resize(nPoints);

    step = (max - min)/nPoints;

    X x = min;

    for(size_t i = 0; i < nPoints; ++i)
    {
        const Y y = Func(x);

        xData[i] = x;
        yData[i] = y;

        x += step;
    }

    PreComputeValues();

    return yData.size();
}


/** AutoSet() : Automatic irregular grid, computed
 *
 *  Irregular grid:
 *    - Pick 25% of the points on a regular grid, for minimum uniform coverage.
 *    - Pick the remaining points from regions of the function with largest second derivative.
 */

template<typename TX, typename TY, int nx, int ny>
size_t  PrecompData<TX, TY, nx, ny>::AutoSet(Y (*Func)(X x), X xmin, X xmax, size_t nPoints)
{
    xMin = xmin;
    xMax = xmax;

    PickBestPoints(Func1, nPoints, overSampling);

    PreComputeValues();

    return yData.size();
}


/// Data retrieval

// Range UNchecked, 0 degree interpolation accessors

template<typename TX, typename TY, int nx, int ny>
typename PrecompData<TX, TY, nx, ny>::Y PrecompData<TX, TY, nx, ny>::operator()(X x) const
{
    return yData[RtoI(x)];
}


// Range checked accessors; check Status()

template<typename TX, typename TY, int nx, int ny>
typename PrecompData<TX, TY, nx, ny>::Y PrecompData<TX, TY, nx, ny>::get(X x) {} //+TODO


// Range checked accessors, interpolated; check Status()

template<typename TX, typename TY, int nx, int ny>
typename PrecompData<TX, TY, nx, ny>::Y PrecompData<TX, TY, nx, ny>::Interpolate(X x)
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


template<typename TX, typename TY, int nx, int ny>
void PrecompData<TX, TY, nx, ny>::Interpolation(int order)
{
    interpolation = order;
}


// Range check

template<typename TX, typename TY, int nx, int ny>
int PrecompData<TX, TY, nx, ny>::RangeCheck(X x)
{
    status = 0;

    if(x < xMin)
        status = wrn_x_less_than_min;
    else if(x > xMax)
        status = wrn_x_more_than_max;

    return status;
}


// Get the whole value set

template<typename TX, typename TY, int nx, int ny>
int PrecompData<TX, TY, nx, ny>::Get(std::vector<X> &_xData , std::vector<Y> &_yData) const
{
    _xData = xData;
    _yData = yData;
    return 0;
}


/// Math functions

template<typename TX, typename TY, int nx, int ny>
TY PrecompData<TX, TY, nx, ny>::FirstDerivative(TX x1, TY y1, TX x2, TY y2) const
{
    /// First derivative (central differences):  d1 = [f(x+1) - f(x)] / [(x+1) - x]
    return  (y2 - y1)/(x2 - x1);
}


template<typename TX, typename TY, int nx, int ny>
TY PrecompData<TX, TY, nx, ny>::SecondDerivative(TX x1, TY y1, TX x2, TY y2, TX x3, TY y3) const
{
    /// Second derivative (central differences):  d2 = [f(x-1) - 2f(x) + f(x+1)] / {[(x+1) - (x-1)]/2}^2
    return  (y1 - 2*y2 + y3)/std::pow(0.5f*(x3 - x1), 2);
}


template<typename TX, typename TY, int nx, int ny>
int PrecompData<TX, TY, nx, ny>::PickBestPoints(Y (*Func)(X x), const size_t nPoints, const float overSampling)
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
