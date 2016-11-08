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
	: interpolation(0), status(0), regularGrid(false), overSampling(2.0f)
{
}

template<typename TX, typename TY, int nx, int ny>
PrecompData<TX, TY, nx, ny>::PrecompData(const std::string _funcName)
	: interpolation(0), status(0), regularGrid(false), overSampling(2.0f), funcName(_funcName)
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
    for(int i = 0; i < nx; ++i)
    {
        kRealInt[i] = yData.size()/(max[i] - min[i]);
        kIntReal[i] = 1/kRealInt[i];
    }

    return 0;
}


// Coordinate <--> index transformations

template<typename TX, typename TY, int nx, int ny>
size_t PrecompData<TX, TY, nx, ny>::RtoI(TX x) const     // real --> integer/index
{
    //+CHECK
    if(nx == 1 && ny == 1)
    	return size_t(kRealInt[0]*(x - min[0]));
    //+TODO - Multidimensional case
}

/* //+?
template<typename TX, typename TY, int nx, int ny>
PrecompData<TX, TY, nx, ny>::TX PrecompData<TX, TY, nx, ny>::ItoR(size_t i) const     // integer/index --> real
{
    //+CHECK
    if(nx == 1 && ny == 1)
        return min[0] + kIntReal[0]*TX(i);
    //+TODO - Multidimensional case
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
    //+CHECK
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

    regularGrid = true;

    return yData.size();
}


template<typename TX, typename TY, int nx, int ny>
size_t  PrecompData<TX, TY, nx, ny>::Set(TY      (*Func)(TX x), 
                                         TX      xmin,
                                         TX      xmax,
                                         size_t  nPoints)
{
    static_assert(nx == 1, "Member function valid for one dimesional independent variable, only.");
    static_assert(ny == 1, "Member function valid for one dimesional dependent variable, only.");

    min[0] = xmin;
    max[0] = xmax;
    
    xData.resize(nPoints);
    yData.resize(nPoints);

    step[0] = (max[0] - min[0])/nPoints;

    X x = min;

    for(size_t i = 0; i < nPoints; ++i)
    {
        const TY y0 = Func(x[0]);

        xData[i] = x;
        yData[i][0] = y0;

        x[0] += step[0];
    }

    PreComputeValues();

    regularGrid = true;

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
	min = xmin;
	max = xmax;

	PickBestPoints(Func, nPoints, overSampling);

    PreComputeValues();

    regularGrid = false;

    return yData.size();
}


/// Data retrieval

// Range UNchecked, 0 degree interpolation accessors

template<typename TX, typename TY, int nx, int ny>
typename PrecompData<TX, TY, nx, ny>::Y PrecompData<TX, TY, nx, ny>::operator()(X x) const
{
    return yData[RtoI(x)];
}


template<typename TX, typename TY, int nx, int ny>
TY PrecompData<TX, TY, nx, ny>::operator()(TX x) const
{
    //+ X xv; xv[0] = x;
    //+ const Y yv = yData[RtoI(xv)];

    const size_t i = RtoI(x);
    const Y yv = yData[i];

    return yv[0];
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

	const TX x0 = ItoR(i);
	const TX x1 = ItoR(i + 1);

	const TY y = yData[i] + (yData[i + 1] - yData[i])*(x - x0)/(x1 - x0);

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

	if(x < min)
        status = wrn_x_less_than_min;
	else if(x > max)
        status = wrn_x_more_than_max;

    return status;
}


template<typename TX, typename TY, int nx, int ny>
int PrecompData<TX, TY, nx, ny>::RangeCheck(TX x)
{
    X xv;
    xv[0] = x;

    return RangeCheck(xv);
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
		TX x; TY y;
		Point(TX _x = 0.0, TY _y = 0.0) : x(_x), y(_y) {}
        bool operator< (const Point &p) const { return x < p.x; }
        bool operator> (const Point &p) const { return x > p.x; }
    };

    struct PointCurv {     // abscissa and second derivative
		TX x; TY d2;
		PointCurv(TX _x = 0.0, TY _d2 = 0.0) : x(_x), d2(_d2) {}
        bool operator> (const PointCurv &p) const { return fabs(d2) > fabs(p.d2); }
    };

    /// Oversample

    const size_t nSamples = size_t(overSampling*nPoints);
	const TX step = (max - min)/nSamples;

    std::vector<PointCurv> samples;
    samples.reserve(nSamples);

	TX x1 = min;
	TX x2 = x1 + step;
	TX x3 = x2 + step;
	TY y1 = Func(x1);
	TY y2 = Func(x2);
	TY y3 = Func(x3);

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

	points[0].x = min;               // always get the first point

    for(size_t i = 1, j = 0; i < nPoints - 1; ++i, ++j)
        points[i].x = samples[j].x;

	points[nPoints - 1].x = max;     // always get the last point

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
