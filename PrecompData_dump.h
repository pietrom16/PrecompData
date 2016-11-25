/**  PrecompData_dump.h

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


/** Dump()
 *  Dump internal values on stdout.
 *  n : if = 0, dump everything;
 *      if > 0, dump the first n points;
 *      if < 0, dump |n| random points.
 */

template<typename TX, typename TY, int nx, int ny>
int PrecompData<TX, TY, nx, ny>::Dump(int n) const
{
	using std::cout;
	using std::endl;

	if(n > int(xData.size()) || n == 0)
		n = xData.size();

	cout << "\n\n----------------------------------------------------------------\n";
	cout << "Dump: " << funcName << "\n" << comment << "\nStatus: " << status << endl;
	cout << "Regular grid: " << regularGrid << "\nOversampling: " << overSampling << endl;

	cout << "\nDim   Min   Max   Step   kRealInt   kIntReal" << endl;
	for(size_t i = 0; i < nx; ++i)
	{
		cout << i + 1 << "   " << min[i] << "   " << max[i] << "   " << step[i] << "   " << kRealInt[i] << "   " << kIntReal[i] << endl;
	}

	if(n > 0)			// dump all or the first n points
	{
		if(n == int(xData.size()))
			cout << "Dump all points:" << endl;
		else
			cout << "Dump the first " << n << " points:" << endl;

		for(size_t j = 0; j < n; ++j)
			DumpElement(j);
	}
	else				// dump |n| random points
	{
		const size_t m = size_t(-n);
		size_t i;

		std::random_device                     rd;
		std::mt19937                           gen(rd());
		std::uniform_int_distribution<size_t>  dist(0, xData.size());

		cout << "Dump " << m << " random points:" << endl;

		for(size_t j = 0; j < m; ++j)
		{
			i = dist(gen);
			DumpElement(i);
		}
	}

	if(xData.size() != yData.size())
	{
		cout << "Warning: size of independent and dependent data do not match." << endl;
		for(size_t j = 0; j < yData.size(); ++j)
		{
			for(size_t i = 0; i < ny; ++i)
			{
				cout << "  " << yData[j][i];
			}
			cout << endl;
		}
	}

	cout << "Dump complete.\n";
	cout << "----------------------------------------------------------------\n";

	return 0;
}


/** DumpElement()
 *  Dump a specific point on stdout.
 *  j : index of the point to dump.
 */

template<typename TX, typename TY, int nx, int ny>
int PrecompData<TX, TY, nx, ny>::DumpElement(size_t j) const
{
	using std::cout;
	using std::endl;

	assert(j < xData.size());

	for(size_t i = 0; i < nx; ++i)
	{
		cout << "  " << xData[j][i];
	}

	if(xData.size() == yData.size())
	{
		cout << " --> ";

		for(size_t i = 0; i < ny; ++i)
		{
			cout << "  " << yData[j][i];
		}
	}

	cout << endl;

	return 0;
}


} // Utilities
