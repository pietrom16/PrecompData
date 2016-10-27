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

float TestFuncLin(float x) {            //   y = 2x
    return 2*x;
}

float TestFuncNonLin1(float x) {        //   y = |x|
    return fabs(x);
}

float TestFuncNonLin2(float x) {        //   y = 1/(|x-2| + 0.1)
    return 1/(fabs(x - 2.0) + 0.1);
}


int main()
{
	using namespace Utilities;

	{ // Test 1 - Interpolation
        cout << "\n\nTest 1: Zero-degree (nearest-neighbor/point sampling/Voronoi) interpolation:" << endl;
		const string funcName = "TestFunc";
		PrecompData<float> itp(funcName);
		itp.SetComment("TestFunc approximation");
		const float x0 = 0.0f, x1 = 6.28f;
		const int nValues = 10;
		const float step = 0.5*(x1 - x0)/nValues;
		itp.Set(&TestFunc, x0, x1, nValues);
		float x = x0;
        float err = 0.0;
		for(int i = 0; i < nValues; ++i) {
            const float y = itp(x);
            err += fabs(TestFunc(x) - y);
			cout << i << ":\t" << funcName << "(" << x << ") = " << TestFunc(x) << " ~ " << y << endl;
			x += step;
		}
        cout << "Total error = " << err << endl;
	}

    { // Test 2 - Interpolation
        cout << "\n\nTest 2: Linear interpolation:" << endl;
        const string funcName = "TestFunc";
        PrecompData<float> itp(funcName);
        const float x0 = 0.0f, x1 = 6.28f;
        const int nValues = 10;
        const float step = 0.5*(x1 - x0)/nValues;
        itp.Set(&TestFunc, x0, x1, nValues);
        float x = x0;
        float err = 0.0;
        itp.Interpolation(1);
        cout << "Interpolation: " << itp.Interpolation() << endl;
        for(int i = 0; i < nValues; ++i) {
            const float y = itp.Interpolate(x);
            err += fabs(TestFunc(x) - y);
            cout << i << ":\t" << funcName << "(" << x << ") = " << TestFunc(x) << " ~ " << y << endl;
            x += step;
        }
        cout << "Total error = " << err << endl;
    }

    { // Test 3 - AutoSet
        cout << "\n\nTest 3: Automatic irregular grid:" << endl;
    }

    return 0;
}

