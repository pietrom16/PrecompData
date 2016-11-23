/**  PrecompData_impl.h

	Copyright 2016 Pietro Mele
	Apache License 2.0
*/

#include "StlExt.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <random>

namespace Utilities {

using namespace stdExt;

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
size_t PrecompData<TX, TY, nx, ny>::VectorToIndex(X x) const      // vector --> index
{
    //+TODO
	/*
		- Thorough search?
		- Mapping index <--> coordinate?
		- Direct computation?
	*/

	//+SLOW - Thorough search

	TX      error     = 0.0f;		// for a specific point
	TX      minError  = std::numeric_limits<TX>::max();
	size_t  minErrPos = -1;

	for(size_t p = 0; p < xData.size(); ++p)
	{
		error = 0.0f;

		for(size_t i = 0; i < nx; ++i)
		{
			TX delta = xData[p][i] - x[i];		//+TEST
			error += delta*delta;
		}

		if(error < minError) {
			minError = error;
			minErrPos = p;
		}
	}

	return minErrPos;
}


template<typename TX, typename TY, int nx, int ny>
size_t PrecompData<TX, TY, nx, ny>::ScalarToIndex(TX x) const     // scalar --> index
{
    static_assert(nx == 1, "Member function valid for one dimesional independent variable, only.");
    return size_t(kRealInt[0]*(x - min[0]));
}


template<typename TX, typename TY, int nx, int ny>
typename PrecompData<TX, TY, nx, ny>::X PrecompData<TX, TY, nx, ny>::IndexToVector(size_t i) const      // index  --> vector
{
    //+TODO
}


template<typename TX, typename TY, int nx, int ny>
TX PrecompData<TX, TY, nx, ny>::IndexToScalar(size_t i) const     // index --> scalar
{
    static_assert(nx == 1, "Member function valid for one dimesional independent variable, only.");
    return min[0] + kIntReal[0]*TX(i);
}



/// Data loading

// Regular grid, computed

template<typename TX, typename TY, int nx, int ny>
size_t  PrecompData<TX, TY, nx, ny>::Set(Y       (*Func)(X x), 
                                         X       xmin,
                                         X       xmax,
                                         size_t  nPoints)
{
    size_t  nSteps[nx];

    // Init
    {
        assert(Func != 0);
        assert(nPoints > 0);
        for(int i = 0; i < nx; ++i)
            assert(xmin[i] < xmax[i]);

        min = xmin;
        max = xmax;

        xData.clear();
        yData.clear();

        FuncX = Func;
        FuncTX = 0;
    }

    // Find step, with these constraints: nPoints, nx, min, max
    // In this context, step is the same across all dimensions.
    {
        TX domainVolume = 1.0f;
        for(int j = 0; j < nx; ++j)
            domainVolume *= (max[j] - min[j]);

        assert(domainVolume > 0.0f);

        const TX elementVolume = domainVolume/nPoints;
    
        TX uniformStep = std::pow(elementVolume, 1/nx);

        for(int j = 0; j < nx; ++j)
        {
            step[j] = uniformStep;
            nSteps[j] = (max[j] - min[j])/step[j];
        }
    }

    //+TEST
    // Scan the nx-dimensional hyperspace; store the computed values
    {
        X x;

        for(size_t i = 0; i < nPoints; ++i)
        {
            size_t dimProd = 1;

            // Transform  i --> X
            for(size_t j = 0; j < nx; ++j)
            {
                x[j] = i / dimProd % nSteps[j];

                dimProd *= nSteps[j];
            }

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


template<typename TX, typename TY, int nx, int ny>
size_t  PrecompData<TX, TY, nx, ny>::Set(TY      (*Func)(TX x), 
                                         TX      xmin,
                                         TX      xmax,
                                         size_t  nPoints)
{
    static_assert(nx == 1, "Member function valid for one dimesional independent variable, only.");
    static_assert(ny == 1, "Member function valid for one dimesional dependent variable, only.");

    FuncTX = Func;
    FuncX = 0;

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
    FuncX = Func;
    FuncTX = 0;

	min = xmin;
	max = xmax;

	PickBestPoints(Func, nPoints, overSampling);

    PreComputeValues();

    regularGrid = false;

    return yData.size();
}


template<typename TX, typename TY, int nx, int ny>
size_t  PrecompData<TX, TY, nx, ny>::AutoSet(TY (*Func)(TX x), TX xmin, TX xmax, size_t nPoints)
{
    static_assert(nx == 1, "Member function valid for one dimesional independent variable, only.");
    static_assert(ny == 1, "Member function valid for one dimesional dependent variable, only.");
    
    FuncTX = Func;
    FuncX = 0;

    min[0] = xmin;
    max[0] = xmax;

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
	return yData[VectorToIndex(x)];
}


template<typename TX, typename TY, int nx, int ny>
TY PrecompData<TX, TY, nx, ny>::operator()(TX x) const
{
    static_assert(nx == 1, "Member function valid for one dimesional independent variable, only.");
    static_assert(ny == 1, "Member function valid for one dimesional dependent variable, only.");

	const size_t i = ScalarToIndex(x);
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


template<typename TX, typename TY, int nx, int ny>
TY PrecompData<TX, TY, nx, ny>::Interpolate(TX x)
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

	//+TODO - Operators < > do not exist for multidimensional data types
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


template<typename TX, typename TY, int nx, int ny>
int PrecompData<TX, TY, nx, ny>::Get(std::vector<TX> &_xData , std::vector<TY> &_yData) const
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


// Evaluate error


// Error on each dimension on known data

template<typename TX, typename TY, int nx, int ny>
typename PrecompData<TX, TY, nx, ny>::Y PrecompData<TX, TY, nx, ny>::EvaluateErrorKnownData() const
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

template<typename TX, typename TY, int nx, int ny>
TY PrecompData<TX, TY, nx, ny>::EvaluateAbsErrorKnownData() const
{
	//+TEST

	const Y error = EvaluateErrorKnownData();
	return Norm(error);
}


// Error on each dimension on random points

template<typename TX, typename TY, int nx, int ny>
typename PrecompData<TX, TY, nx, ny>::Y PrecompData<TX, TY, nx, ny>::EvaluateError(int nTestPoints) const
{
	//+TEST

	X x;
	Y error;

	std::random_device  rd;
	std::minstd_rand    gen(rd());

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

template<typename TX, typename TY, int nx, int ny>
TY PrecompData<TX, TY, nx, ny>::EvaluateAbsError(int nTestPoints) const
{
	//+TEST

	const Y error = EvaluateError(nTestPoints);
	return Norm(error);
}


/// Math functions


template<typename TX, typename TY, int nx, int ny>
TY PrecompData<TX, TY, nx, ny>::Norm(const Y &y) const
{
    TY norm = 0.0;
    
    for(size_t i = 0; i < ny; ++i)
        norm += y[i] * y[i];

    norm = std::sqrt(norm);

    return norm;
}


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


template<typename TX, typename TY, int nx, int ny>
int PrecompData<TX, TY, nx, ny>::PickBestPoints(TY (*Func)(TX x), const size_t nPoints, const float overSampling)
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
