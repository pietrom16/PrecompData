/**  PrecompData_impl.h

	Copyright 2016 Pietro Mele
	Apache License 2.0
*/

#include "StlExt.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>  //+T+++
#include <limits>
#include <random>

namespace Utilities {

using namespace stdExt;

template<int nPoints, typename TX, typename TY, int ny>
PrecompData<nPoints, TX, TY, ny>::PrecompData()
    : PrecompData_base()
{
}

template<int nPoints, typename TX, typename TY, int ny>
PrecompData<nPoints, TX, TY, ny>::PrecompData(const std::string _funcName)
    : PrecompData_base(_funcName)
{
}


// Precompute constant values

template<int nPoints, typename TX, typename TY, int ny>
int PrecompData<nPoints, TX, TY, ny>::PreComputeValues()
{
    //+CHECK
    // Set up conversion constants
	kRealInt = yData.size()/(max - min);
	kIntReal = 1/kRealInt;

    return 0;
}


// Coordinate <--> index transformations


template<int nPoints, typename TX, typename TY, int ny>
size_t PrecompData<nPoints, TX, TY, ny>::ScalarToIndex(TX x) const     // scalar --> index
{
	return size_t(kRealInt*(x - min));
}


template<int nPoints, typename TX, typename TY, int ny>
TX PrecompData<nPoints, TX, TY, ny>::IndexToScalar(size_t i) const     // index --> scalar
{
	return min + kIntReal*TX(i);
}



/// Data loading

// Regular grid, computed

template<int nPoints, typename TX, typename TY, int ny>
size_t  PrecompData<nPoints, TX, TY, ny>::Set(YData  (*Func)(TX x),
                                              TX      xmin,
                                              TX      xmax)
{
    // Init
    {
        assert(Func != 0);
        assert(nPoints > 0);
		assert(xmin < xmax);

        min = xmin;
        max = xmax;

        xData.clear();
        yData.clear();
		xData.reserve(nPoints);
		yData.reserve(nPoints);

		FuncTXVY = Func;
		FuncTXTY = 0;
    }

	// Find step, with these constraints: nPoints, min, max
    {
		const TX domainInterval = (max - min);

		assert(domainInterval > 0.0f);

		step = domainInterval/nPoints;
    }

    //+TEST
	// Scan the interval on the x axis; store the computed values
    {
        for(size_t i = 0; i < nPoints; ++i)
        {
			// Transform  i --> x
			const TX x = step*i + xmin;

            const Y y = Func(x);

            xData.push_back(x);
            yData.push_back(y);

            // Check independent and dependent vectors are aligned
            assert(xData.size() == yData.size());
        }
    }

    PreComputeValues();

    regularGrid = true;

    return yData.size();
}


template<int nPoints, typename TX, typename TY, int ny>
size_t  PrecompData<nPoints, TX, TY, ny>::Set(TY  (*Func)(TX x),
                                              TX  xmin,
                                              TX  xmax)
{
    static_assert(ny == 1, "Member function valid for one dimesional dependent variable, only.");

	// Init
	{
		assert(Func != 0);
		assert(nPoints > 0);
		assert(xmin < xmax);

		min = xmin;
		max = xmax;

		xData.resize(nPoints);
		yData.resize(nPoints);

		FuncTXVY = 0;
		FuncTXTY = Func;
	}

	// Find step, with these constraints: nPoints, min, max
	{
		const TX domainInterval = (max - min);

		assert(domainInterval > 0.0f);

		step = domainInterval/nPoints;
	}

	//+TEST
	// Scan the interval on the x axis; store the computed values
	{
		for(size_t i = 0; i < nPoints; ++i)
		{
			// Transform  i --> x
			const TX x = step*i + xmin;

			const TY y = Func(x);

			xData.push_back(x);
			yData.push_back(y);

			// Check independent and dependent vectors are aligned
			assert(xData.size() == yData.size());
		}
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

template<int nPoints, typename TX, typename TY, int ny>
size_t  PrecompData<nPoints, TX, TY, ny>::AutoSet(YData (*Func)(TX x), TX xmin, TX xmax)
{
	FuncTXVY = Func;
	FuncTXTY = 0;

	min = xmin;
	max = xmax;

	PickBestPoints(Func, nPoints, overSampling);

    PreComputeValues();

    regularGrid = false;

    return yData.size();
}


template<int nPoints, typename TX, typename TY, int ny>
size_t  PrecompData<nPoints, TX, TY, ny>::AutoSet(TY (*Func)(TX x), TX xmin, TX xmax)
{
    static_assert(ny == 1, "Member function valid for one dimesional dependent variable, only.");
    
	FuncTXVY = 0;
	FuncTXTY = Func;

	min = xmin;
	max = xmax;

    PickBestPoints(Func, nPoints, overSampling);

    PreComputeValues();

    regularGrid = false;

    return yData.size();
}


/// Data retrieval

// Range UNchecked, 0 degree interpolation accessors

template<int nPoints, typename TX, typename TY, int ny>
size_t PrecompData<nPoints, TX, TY, ny>::operator()(TX _x, TY &_y) const
{
	static_assert(ny == 1, "Member function valid for one dimesional dependent variable, only.");

	const size_t i = ScalarToIndex(_x);
	_y = yData[i];

	return i;
}


template<int nPoints, typename TX, typename TY, int ny>
size_t PrecompData<nPoints, TX, TY, ny>::operator()(TX _x, YData &_y) const
{
	const size_t i = VectorToIndex(_x);
	_y = yData[i];

	return i;
}


template<int nPoints, typename TX, typename TY, int ny>
TY  PrecompData<nPoints, TX, TY, ny>::operator()(TX _x) const
{
	static_assert(ny == 1, "Member function valid for one dimesional dependent variable, only.");

	const size_t i = ScalarToIndex(_x);
	return yData[i];
}


// Range checked, 0 degree interpolation accessors

template<int nPoints, typename TX, typename TY, int ny>
size_t PrecompData<nPoints, TX, TY, ny>::get(TX _x, TY &_y) const
{
	if(_x < min) return wrn_x_less_than_min;
	if(_x > max) return wrn_x_more_than_max;

	return operator(_x, _y);
}

template<int nPoints, typename TX, typename TY, int ny>
size_t PrecompData<nPoints, TX, TY, ny>::get(TX _x, YData &_y) const
{
	if(_x < min) return wrn_x_less_than_min;
	if(_x > max) return wrn_x_more_than_max;

	return operator(_x, _y);
}


// Range checked accessors, interpolated; check Status()

template<int nPoints, typename TX, typename TY, int ny>
typename PrecompData<nPoints, TX, TY, ny>::Y PrecompData<nPoints, TX, TY, ny>::Interpolate(X x)
{
    RangeCheck(x);

	const size_t i = VectorToIndex(x);

	if(i >= xData.size() - 2) {
		status = wrn_x_more_than_max;
        return yData.back();
	}

	const X x0 = IndexToVector(i);
	const X x1 = IndexToVector(i + 1);

    //+TODO: Fix - multidimensional
	const TY y = yData[i] + (yData[i + 1] - yData[i])*(x - x0)/(x1 - x0);

    return y;
}


template<int nPoints, typename TX, typename TY, int ny>
TY PrecompData<nPoints, TX, TY, ny>::Interpolate(TX x)
{
    static_assert(nx == 1, "Member function valid for one dimesional independent variable, only.");
    static_assert(ny == 1, "Member function valid for one dimesional dependent variable, only.");

    RangeCheck(x);

	const size_t i = ScalarToIndex(x);

	if(i >= xData.size() - 2) {
		status = wrn_x_more_than_max;
		return yData.back()[0];
	}

	const TX x0 = IndexToScalar(i);
	const TX x1 = IndexToScalar(i + 1);

    const TY yv = yData[i][0] + (yData[i + 1][0] - yData[i][0])*(x - x0)/(x1 - x0);

    return yv;
}


// Range check

template<int nPoints, typename TX, typename TY, int ny>
int PrecompData<nPoints, TX, TY, ny>::RangeCheck(TX x)
{
	if(x < min) return -1;
	if(x > max) return 1;

	return 0;	// OK, in range
}


// Get the whole value set

template<int nPoints, typename TX, typename TY, int ny>
int PrecompData<nPoints, TX, TY, ny>::Get(std::vector<X> &_xData , std::vector<Y> &_yData) const
{
    _xData = xData;
    _yData = yData;
    return 0;
}


template<int nPoints, typename TX, typename TY, int ny>
int PrecompData<nPoints, TX, TY, ny>::Get(std::vector<TX> &_xData , std::vector<TY> &_yData) const
{
    static_assert(nx == 1, "Member function valid for one dimesional independent variable, only.");
    static_assert(ny == 1, "Member function valid for one dimesional dependent variable, only.");
    
    //+TEST

    _xData.clear();
    
    for(size_t i = 0; i < xData.size(); ++i)
    {
        TX p = xData[i][0];     //+? Swap indices?
        _xData.push_back(p);
    }

    _yData.clear();

    for(size_t i = 0; i < yData.size(); ++i)
    {
        TY p = yData[i][0];     //+? Swap indices?
        _yData.push_back(p);
    }

//+B    _xData = xData[0];
//+B    _yData = yData[0];

    return 0;
}


/// Evaluate error


// Error on each dimension on known data

template<int nPoints, typename TX, typename TY, int ny>
typename PrecompData<nPoints, TX, TY, ny>::Y PrecompData<nPoints, TX, TY, ny>::EvaluateErrorKnownData() const
{
	//+TEST

	Y error;
    
    for(size_t j = 0; j < error.size(); ++j)
        error[j] = 0.0;

    for(size_t i = 0; i < xData.size(); ++i)
    {
        const X x      = xData[i];
        const Y y      = yData[i];
        const Y y_comp = FuncX(x);

        for(size_t j = 0; j < error.size(); ++j) {
            //error[j] += fabs(y_comp[j] - y[j]);                  // mean absolute error
            error[j] += (y_comp[j] - y[j])*(y_comp[j] - y[j]);     // mean squared error
        }
    }

    for(size_t j = 0; j < error.size(); ++j)
        error[j] /= xData.size();

    return error;
}


// Absolute error on known data

template<int nPoints, typename TX, typename TY, int ny>
TY PrecompData<nPoints, TX, TY, ny>::EvaluateAbsErrorKnownData() const
{
	//+TEST

	const Y error = EvaluateErrorKnownData();
	return Norm(error);
}


// Error on each dimension on random points

template<int nPoints, typename TX, typename TY, int ny>
typename PrecompData<nPoints, TX, TY, ny>::Y PrecompData<nPoints, TX, TY, ny>::EvaluateError(int nTestPoints) const
{
	//+TEST

	X x;
	Y error;

	std::random_device  rd;
	std::mt19937        gen(rd());

	for(size_t j = 0; j < error.size(); ++j)
		error[j] = 0.0;

	for(size_t i = 0; i < nTestPoints; ++i)
	{
		// Set a random x
		for(size_t j = 0; j < x.size(); ++j)
		{
			std::uniform_real_distribution<> dist(min[j], max[j]);
			x[j] = dist(gen);
		}

		const Y y      = Interpolate(x);
		const Y y_comp = FuncX(x);

		for(size_t j = 0; j < error.size(); ++j) {
			//error[j] += fabs(y_comp[j] - y[j]);                  // mean absolute error
			error[j] += (y_comp[j] - y[j])*(y_comp[j] - y[j]);     // mean squared error
		}
	}

	for(size_t j = 0; j < error.size(); ++j)
		error[j] /= xData.size();

	return error;
}


// Absolute error on random points

template<int nPoints, typename TX, typename TY, int ny>
TY PrecompData<nPoints, TX, TY, ny>::EvaluateAbsError(int nTestPoints) const
{
	//+TEST

	const Y error = EvaluateError(nTestPoints);
	return Norm(error);
}


/// Math functions


template<int nPoints, typename TX, typename TY, int ny>
TY PrecompData<nPoints, TX, TY, ny>::Norm(const Y &y) const
{
    TY norm = 0.0;
    
    for(size_t i = 0; i < ny; ++i)
        norm += y[i] * y[i];

    norm = std::sqrt(norm);

    return norm;
}


template<int nPoints, typename TX, typename TY, int ny>
TY PrecompData<nPoints, TX, TY, ny>::FirstDerivative(TX x1, TY y1, TX x2, TY y2) const
{
    /// First derivative (central differences):  d1 = [f(x+1) - f(x)] / [(x+1) - x]
    return  (y2 - y1)/(x2 - x1);
}


template<int nPoints, typename TX, typename TY, int ny>
TY PrecompData<nPoints, TX, TY, ny>::SecondDerivative(TX x1, TY y1, TX x2, TY y2, TX x3, TY y3) const
{
    /// Second derivative (central differences):  d2 = [f(x-1) - 2f(x) + f(x+1)] / {[(x+1) - (x-1)]/2}^2
    return  (y1 - 2*y2 + y3)/std::pow(0.5f*(x3 - x1), 2);
}


template<int nPoints, typename TX, typename TY, int ny>
int PrecompData<nPoints, TX, TY, ny>::PickBestPoints(Y (*Func)(X x), const float overSampling)
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


template<int nPoints, typename TX, typename TY, int ny>
int PrecompData<nPoints, TX, TY, ny>::PickBestPoints(TY (*Func)(TX x), const float overSampling)
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
    const TX step = (max[0] - min[0])/nSamples;

    std::vector<PointCurv> samples;
    samples.reserve(nSamples);

    TX x1 = min[0];
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
        y3 = Func(x3);
    }

    // Sort based on decreasing second derivative absolute value
    std::sort(samples.begin(), samples.end(), std::greater<PointCurv>());

    /// Pick the points with highest second derivative (curvature)

    std::vector<Point> points;
    points.resize(nPoints);

    points[0].x = min[0];               // always get the first point

    for(size_t i = 1, j = 0; i < nPoints - 1; ++i, ++j)
        points[i].x = samples[j].x;

    points[nPoints - 1].x = max[0];     // always get the last point

    for(size_t i = 0; i < nPoints; ++i)
        points[i].y = Func(points[i].x);

    // Sort the picked points based on increasing abscissa
    std::sort(points.begin(), points.end());

    // Copy data in member variables

    xData.clear();
    yData.clear();
    xData.resize(points.size());
    yData.resize(points.size());

    for(size_t i = 0; i < points.size(); ++i)
    {
        xData[i][0] = points[i].x;
        yData[i][0] = points[i].y;
    }

    return 0;
}


} // Utilities
