/**  StlExt.hpp

	Extensions to the Standard Template Library:
	std::array: operator- and operator/

	Copyright 2016 Pietro Mele
	Apache License 2.0
*/

#ifndef STLEXT_HPP
#define STLEXT_HPP

#include <array>
#include <cassert>

namespace stdExt {


/// --- std::array ---


template<class T, std::size_t N>
std::array<T,N> operator-(const std::array<T,N>& lhs,
                          const std::array<T,N>& rhs)
{
	std::array<T,N> r;
	for(size_t i = 0; i < N; ++i)
		r[i] = lhs[i] - rhs[i];
	return r;
}


template<class T, std::size_t N>
std::array<T,N> operator/(const std::array<T,N>& lhs,
                          const T rhs)
{
	std::array<T,N> r;
	for(size_t i = 0; i < N; ++i)
		r[i] = lhs[i]/rhs;
	return r;
}


} // stdExt

#endif // STLEXT_HPP
