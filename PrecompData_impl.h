/**  PrecompData_impl.h

	Copyright 2016 Pietro Mele
	Apache License 2.0
*/

#include "PrecompData.h"

namespace Utilities {

template<typename T>
PrecompData<T>::PrecompData()
	: interpolation(0), status(0)
{
}

template<typename T>
PrecompData<T>::PrecompData(const std::string _funcName)
	: interpolation(0), status(0), funcName(_funcName)
{
}


template<typename T>
int  PrecompData<T>::SetFunctionName(const std::string &_funcName)
{
	funcName = _funcName;
	return 0;
}

template<typename T>
int  PrecompData<T>::SetComment(const std::string &_comment)
{
	comment = _comment;
	return 0;
}

template<typename T>
std::string  PrecompData<T>::FunctionName() const
{
	return funcName;
}

template<typename T>
std::string  PrecompData<T>::Comment() const
{
	return comment;
}


// Regular grid

template<typename T>
int  PrecompData<T>::Set(T (*Func1)(T x), T xmin, T xmax, size_t nPoints)     // line
{
	line.resize(nPoints);
	const T step = (xmax - xmin)/nPoints;

	T x = xmin;

	for(size_t i = 0; i < nPoints; ++i)
	{
		const T y = Func1(x);

		line[i] = y;

		x += step;
	}

	return line.size();
}


template<typename T>
int  PrecompData<T>::Set(T (*Func2)(T x, T y),
                      T xmin, T xmax, size_t xnPoints,
                      T ymin, T ymax, size_t ynPoints);    // plane

template<typename T>
int  PrecompData<T>::Set(T (*Func3)(T x, T y, T z),
                      T xmin, T xmax, size_t xnPoints,
                      T ymin, T ymax, size_t ynPoints,
                      T zmin, T zmax, size_t znPoints);    // volume

// Irregular grid

template<typename T>
int  PrecompData<T>::Set(T (*Func1)(T x), const std::vector<T> &x);      // line

template<typename T>
int  PrecompData<T>::Set(T (*Func2)(T x, T y),
                      const std::vector<T> &x,
                      const std::vector<T> &y);      // plane

template<typename T>
int  PrecompData<T>::Set(T (*Func3)(T x, T y, T z),
                      const std::vector<T> &x,
                      const std::vector<T> &y,
                      const std::vector<T> &z);      // volume

// Range UNchecked accessors

template<typename T>
T PrecompData<T>::operator()(T x) const;

template<typename T>
T PrecompData<T>::operator()(T x, T y) const;

template<typename T>
T PrecompData<T>::operator()(T x, T y, T z) const;

// Range checked accessors; check Status()

template<typename T>
T PrecompData<T>::get(T x) const;

template<typename T>
T PrecompData<T>::get(T x, T y) const;

template<typename T>
T PrecompData<T>::get(T x, T y, T z) const;

// Range checked accessors, interpolated; check Status()

template<typename T>
T PrecompData<T>::Interpolate(T x) const;

template<typename T>
T PrecompData<T>::Interpolate(T x, T y) const;

template<typename T>
T PrecompData<T>::Interpolate(T x, T y, T z) const;


template<typename T>
void PrecompData<T>::Interpolate(int order);


template<typename T>
int PrecompData<T>::CopyOnGPU();

// Copy a subset

template<typename T>
int PrecompData<T>::CopyOnGPU(T xbeg, T xend);

template<typename T>
int PrecompData<T>::CopyOnGPU(T xbeg, T xend, T ybeg, T yend);

template<typename T>
int PrecompData<T>::CopyOnGPU(T xbeg, T xend, T ybeg, T yend, T zbeg, T zend);

} // Utilities
