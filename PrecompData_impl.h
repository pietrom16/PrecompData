/**  PrecompData_impl.h

	Copyright 2016 Pietro Mele
	Apache License 2.0
*/

#include "StlExt.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <chrono>
#include <functional>
#include <iostream>  //+T+++
#include <limits>
#include <random>
#include <vector>

namespace Utilities {

using namespace stdExt;

template<int nPoints, typename TX, typename TY>
PrecompData<nPoints, TX, TY>::PrecompData()
    : PrecompData_base()
{
}

template<int nPoints, typename TX, typename TY>
PrecompData<nPoints, TX, TY>::PrecompData(const std::string _funcName)
    : PrecompData_base(_funcName)
{
}


// Precompute constant values

template<int nPoints, typename TX, typename TY>
int PrecompData<nPoints, TX, TY>::PreComputeValues()
{
    //+CHECK
    // Set up conversion constants
	kRealInt = xData.size()/(xMax - xMin);
	kIntReal = 1/kRealInt;

    return 0;
}


// Coordinate <--> index transformations


template<int nPoints, typename TX, typename TY>
size_t PrecompData<nPoints, TX, TY>::ScalarToIndex(TX x) const     // scalar --> index
{
	return size_t(kRealInt*(x - xMin));
}


template<int nPoints, typename TX, typename TY>
TX PrecompData<nPoints, TX, TY>::IndexToScalar(size_t i) const     // index --> scalar
{
	return xMin + kIntReal*TX(i);
}



/// Data loading

// Regular grid, computed

template<int nPoints, typename TX, typename TY>
size_t  PrecompData<nPoints, TX, TY>::set(TY  (*Func)(TX x),
                                          TX  xmin,
                                          TX  xmax)
{
	// Init
	{
		assert(Func != 0);
		assert(nPoints > 0);
		assert(xmin < xmax);

		xMin = xmin;
		xMax = xmax;

		FuncXY = Func;
	}

	// Find step, with these constraints: nPoints, min, max
	{
		const TX domainInterval = (xMax - xMin);

		assert(domainInterval > 0.0f);

		step = domainInterval/nPoints;
	}

	//+TEST
	// Scan the interval on the x axis; store the computed values
	{
		// Check independent and dependent vectors are aligned
		assert(xData.size() == yData.size());

		for(size_t i = 0; i < nPoints; ++i)
		{
			// Transform  i --> x
			const TX x = step*i + xMin;

			const TY y = FuncXY(x);

			xData[i] = x;
			yData[i] = y;
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

template<int nPoints, typename TX, typename TY>
size_t  PrecompData<nPoints, TX, TY>::AutoSet(TY (*Func)(TX x), TX xmin, TX xmax)
{
	FuncXY = Func;

	xMin = xmin;
	xMax = xmax;

	PickBestPoints(FuncXY, overSampling);

    PreComputeValues();

    regularGrid = false;

    return yData.size();
}


/// Data retrieval

// Range UNchecked, 0 degree interpolation accessors

template<int nPoints, typename TX, typename TY>
size_t PrecompData<nPoints, TX, TY>::operator()(TX _x, TY &_y) const
{
	const size_t i = ScalarToIndex(_x);
	_y = yData[i];

	return i;
}


template<int nPoints, typename TX, typename TY>
TY  PrecompData<nPoints, TX, TY>::operator()(TX _x) const
{
	const size_t i = ScalarToIndex(_x);
	return yData[i];
}


// Range checked, 0 degree interpolation accessors

template<int nPoints, typename TX, typename TY>
size_t PrecompData<nPoints, TX, TY>::get(TX _x, TY &_y) const
{
	size_t i = 0;
	if(_x < xMin) { _x = xMin; i = wrn_x_less_than_min; }
	if(_x > xMax) { _x = xMax; i = wrn_x_more_than_max; }

	operator()(_x, _y);

	return i;
}


// Range checked accessors, linear interpolation

template<int nPoints, typename TX, typename TY>
size_t PrecompData<nPoints, TX, TY>::Interpolate(TX _x, TY &_y) const
{
	if(_x <= xMin) {
		_y = yData[0];
		return wrn_x_less_than_min;
	}

	if(_x >= xMax) {
		_y = yData.back();
		return wrn_x_more_than_max;
	}

	const size_t i = ScalarToIndex(_x);

	const TX x0 = IndexToScalar(i);
	const TX x1 = IndexToScalar(i + 1);

	assert(_x >= x0);
	assert(_x <= x1);

	if(i < yData.size())
		_y = yData[i] + (yData[i + 1] - yData[i])*(_x - x0)/(x1 - x0);	//+TEST
	else
		_y = yData.back();

	return i;
}


template<int nPoints, typename TX, typename TY>
TY PrecompData<nPoints, TX, TY>::Interpolate(TX _x)
{
	TY y;

	Interpolate(_x, y);

	return y;
}


// Range check

template<int nPoints, typename TX, typename TY>
int PrecompData<nPoints, TX, TY>::RangeCheck(TX x)
{
	if(x < xMin) return -1;
	if(x > xMax) return 1;

	return 0;	// OK, in range
}


// Get the whole value set

template<int nPoints, typename TX, typename TY>
int PrecompData<nPoints, TX, TY>::get(std::vector<TX> &_xData , std::vector<TY> &_yData) const
{
	if(_xData.capacity() < xData.size())
		_xData.reserve(xData.size());

	if(_yData.capacity() < yData.size())
		_yData.reserve(yData.size());

	std::copy_n(xData.cbegin(), xData.size(), std::back_inserter(_xData));
	std::copy_n(yData.cbegin(), yData.size(), std::back_inserter(_yData));

	return 0;
}


/// Evaluate error


// Error on known data

template<int nPoints, typename TX, typename TY>
TY PrecompData<nPoints, TX, TY>::EvaluateErrorKnownData() const
{
	//+TEST

	assert(FuncXY);
	if(xData.size() == 0) return 0.0;

	TY  error = 0.0;
    
	for(size_t i = 0; i < xData.size(); ++i)
    {
		const TX  x = xData[i];
		const TY  y = yData[i];
		const TY  y_comp = FuncXY(x);

		//error += fabs(y_comp - y);            // mean absolute error
		error += (y_comp - y)*(y_comp - y);     // mean squared error
    }

	error /= xData.size();

    return error;
}


// Absolute error on known data

template<int nPoints, typename TX, typename TY>
TY PrecompData<nPoints, TX, TY>::EvaluateAbsErrorKnownData() const
{
	return Norm(EvaluateErrorKnownData());
}


// Error on random points

template<int nPoints, typename TX, typename TY>
TY PrecompData<nPoints, TX, TY>::EvaluateError(int nTestPoints) const
{
	//+TEST

	assert(FuncXY);
	if(xData.size() == 0) return 0.0;

	TX  x;
	TY  y, error = 0.0;

	std::random_device                rd;
	std::mt19937                      gen(rd());
	std::uniform_real_distribution<>  dist(xMin, xMax);

	for(size_t i = 0; i < nTestPoints; ++i)
	{
		// Set a random x
		x = dist(gen);

		Interpolate(x, y);
		const TY y_comp = FuncXY(x);

		//error += fabs(y_comp - y);            // mean absolute error
		error += (y_comp - y)*(y_comp - y);     // mean squared error
	}

	error /= xData.size();

	return error;
}


// Absolute error on random points

template<int nPoints, typename TX, typename TY>
TY PrecompData<nPoints, TX, TY>::EvaluateAbsError(int nTestPoints) const
{
	return Norm(EvaluateError(nTestPoints));
}


// Performance evaluation

template<int nPoints, typename TX, typename TY>
float PrecompData<nPoints, TX, TY>::PerformanceImprovement(int _nTestPoints)
{
	using namespace std::chrono;

	// Better if > 1, worse if in [0, 1], error if < 0
	float ratio = 0.0;

	if(FuncXY == 0) return -1.0;

	if(_nTestPoints == 0) _nTestPoints = 10*nPoints;

	const TX step = (xMax - xMin)/_nTestPoints;

	duration<float> timeComp;
	duration<float> timePrecomp;

	std::vector<float> results;  // store them to avoid loop optimizations
	results.resize(_nTestPoints);

	// Find time to do real-time computations
	{
		time_point<system_clock> start = system_clock::now();

		TX x = xMin;
		for(int i = 0; i < _nTestPoints; ++i)
		{
			results[i] = FuncXY(x);
			x += step;
		}

		time_point<system_clock> end = system_clock::now();
		timeComp = end - start;
	}

	// Find time to interpolate with precomputations
	{
		time_point<system_clock> start = system_clock::now();

		TX x = xMin;
		for(int i = 0; i < _nTestPoints; ++i)
		{
			results[i] = Interpolate(x);
			x += step;
		}

		time_point<system_clock> end = system_clock::now();
		timePrecomp = end - start;
	}

	std::cout << "timeComp = " << timeComp.count() << "   timePrecomp = " << timePrecomp.count() << std::endl; //+T+++

	ratio = timeComp/timePrecomp;

	return ratio;
}


/// Math functions


template<int nPoints, typename TX, typename TY>
TY PrecompData<nPoints, TX, TY>::Norm(const TY &y) const
{
	return std::fabs(y);
}


template<int nPoints, typename TX, typename TY>
TY PrecompData<nPoints, TX, TY>::FirstDerivative(TX x1, TY y1, TX x2, TY y2) const
{
    /// First derivative (central differences):  d1 = [f(x+1) - f(x)] / [(x+1) - x]
    return  (y2 - y1)/(x2 - x1);
}


template<int nPoints, typename TX, typename TY>
TY PrecompData<nPoints, TX, TY>::SecondDerivative(TX x1, TY y1, TX x2, TY y2, TX x3, TY y3) const
{
    /// Second derivative (central differences):  d2 = [f(x-1) - 2f(x) + f(x+1)] / {[(x+1) - (x-1)]/2}^2
    return  (y1 - 2*y2 + y3)/std::pow(0.5f*(x3 - x1), 2);
}


//+D+
/*
template<int nPoints, typename TX, typename TY>
int PrecompData<nPoints, TX, TY>::PickBestPoints(YData (*Func)(TX x), const float overSampling)
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
	YData y1 = Func(x1);
	YData y2 = Func(x2);
	YData y3 = Func(x3);

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
*/


template<int nPoints, typename TX, typename TY>
int PrecompData<nPoints, TX, TY>::PickBestPoints(TY (*Func)(TX x), const float overSampling)
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
	const TX step = (xMax - xMin)/nSamples;

    std::vector<PointCurv> samples;
    samples.reserve(nSamples);

	TX x1 = xMin;
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

	points[0].x = xMin;               // always get the first point

    for(size_t i = 1, j = 0; i < nPoints - 1; ++i, ++j)
        points[i].x = samples[j].x;

	points[nPoints - 1].x = xMax;     // always get the last point

    for(size_t i = 0; i < nPoints; ++i)
        points[i].y = Func(points[i].x);

    // Sort the picked points based on increasing abscissa
    std::sort(points.begin(), points.end());

    // Copy data in member variables

    for(size_t i = 0; i < points.size(); ++i)
    {
		xData[i] = points[i].x;
		yData[i] = points[i].y;
    }

    return 0;
}


} // Utilities
