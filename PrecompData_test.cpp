/// PrecompData_test.cpp

/** Test the PrecompData class
 */

#include "PrecompData.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

float TestFunc(float x) {
	return sin(x);
}


int main()
{
	using namespace Utilities;

	{ // Test 1
        cout << "\n\nTest 1: Zero-degree (nearest-neighbor/point sampling/Voronoi) interpolation:" << endl;
		const string funcName = "TestFunc";
		PrecompData<float> itp(funcName);
		itp.SetComment("TestFunc approximation");
		const float x0 = 0.0f, x1 = 6.28f;
		const int nValues = 10;
		const float step = 0.5*(x1 - x0)/nValues;
		itp.Set(&TestFunc, x0, x1, nValues);
		float x = x0;
		for(int i = 0; i < nValues; ++i) {
			cout << i << ":\t" << funcName << "(" << x << ") = " << TestFunc(x) << " ~ " << itp(x) << endl;
			x += step;
		}
	}

    { // Test 2
        cout << "\n\nTest 2: Linear interpolation:" << endl;
        const string funcName = "TestFunc";
        PrecompData<float> itp(funcName);
        const float x0 = 0.0f, x1 = 6.28f;
        const int nValues = 10;
        const float step = 0.5*(x1 - x0)/nValues;
        itp.Set(&TestFunc, x0, x1, nValues);
        float x = x0;
        itp.Interpolation(1);
        cout << "Interpolation: " << itp.Interpolation() << endl;
        for(int i = 0; i < nValues; ++i) {
            cout << i << ":\t" << funcName << "(" << x << ") = " << TestFunc(x) << " ~ " << itp.Interpolate(x) << endl;
            x += step;
        }
    }

    return 0;
}

