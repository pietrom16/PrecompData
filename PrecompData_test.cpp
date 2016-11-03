/// PrecompData_test.cpp

/** Test the PrecompData class
 */

#include "PrecompData.h"
#include "PrecompData_test.h"
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
    return 1/(fabs(x - 2.0f) + 0.1f);
}

float TestFuncNonLinSin(float x) {        //   y = sin(x)
    return sin(x);
}


namespace Utilities {

PrecompData_test::PrecompData_test()
{
	using namespace Utilities;

    cout << std::fixed;

    const int nValues = 20;

    // Test 1 - Interpolation
	{
        cout << "\n\nTest 1: Zero-degree (nearest-neighbor/point sampling/Voronoi) interpolation:" << endl;
		const string funcName = "TestFunc";
		PrecompData<float> itp(funcName);
		itp.SetComment("TestFunc approximation");
		const float x0 = 0.0f, x1 = 6.28f;
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

    // Test 2 - Interpolation
    {
        cout << "\n\nTest 2: Linear interpolation:" << endl;
        const string funcName = "TestFunc";
        PrecompData<float> itp(funcName);
        const float x0 = 0.0f, x1 = 6.28f;
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

    // Test 3 - AutoSet:  y = 2x
    {
        cout << "\n\nTest 3: Automatic irregular grid:    y = 2x" << endl;
        const string funcName = "y = 2x";
        PrecompData<float> itp(funcName);
        const float x0 = 0.0f, x1 = 6.28f;
        itp.AutoSet(&TestFuncLin, x0, x1, nValues);
        cerr << "x0 = " << x0 << "  x1 = " << x1 << "  nValues = " << nValues << endl;  //+T+
        std::vector<float> vx, vy;
        itp.Get(vx, vy);
        cout << "Sizes:  x = " << vx.size() << ";  y = " << vy.size() << endl;
        for(size_t i = 0; i < nValues; ++i) {
            cout << i << ":  " << vx[i] << ", " << vy[i] << endl;
        }
    }

    // Test 4 - AutoSet:  y = 1/(|x-2| + 0.1)
    {
        cout << "\n\nTest 4: Automatic irregular grid:    y = 1/(|x-2| + 0.1)" << endl;
        const string funcName = "y = 1/(|x-2| + 0.1)";
        PrecompData<float> itp(funcName);
        const float x0 = 0.0f, x1 = 6.28f;
        itp.AutoSet(&TestFuncNonLin2, x0, x1, nValues);
        cerr << "x0 = " << x0 << "  x1 = " << x1 << "  nValues = " << nValues << endl;  //+T+
        std::vector<float> vx, vy;
        itp.Get(vx, vy);
        cout << "Sizes:  x = " << vx.size() << ";  y = " << vy.size() << endl;
        for(size_t i = 0; i < nValues; ++i) {
            cout << i << ":  " << vx[i] << ", " << vy[i] << endl;
        }
    }

    // Test 5 - Derivatives
    {
        cout << "\n\nTest 5: Derivatives" << endl;
        int nTests = 0, nFailed = 0;
        const string funcName = "Derivatives";
        float x1, y1, x2, y2, x3, y3, der1, der2, expRes;
        PrecompData<float> test;

        // First derivative
        x1 = 0.0; y1 = 0.0; x2 = 1.0; y2 = 0.0; expRes = 0.0;
        der1 = test.FirstDerivative(x1, y1, x2, y2);
        ++nTests;
        if(fabs(der1 - expRes) > 0.0001) {
            ++nFailed;
            cerr << "Error - First derivative 1" << endl;
        }
        x1 = 0.0; y1 = 0.0; x2 = 1.0; y2 = 1.0; expRes = 1.0;
        der1 = test.FirstDerivative(x1, y1, x2, y2);
        ++nTests;
        if(fabs(der1 - expRes) > 0.0001) {
            ++nFailed;
            cerr << "Error - First derivative 2" << endl;
        }
        x1 = 1.0; y1 = 0.0; x2 = 0.0; y2 = 1.0; expRes = -1.0;
        der1 = test.FirstDerivative(x1, y1, x2, y2);
        ++nTests;
        if(fabs(der1 - expRes) > 0.0001) {
            ++nFailed;
            cerr << "Error - First derivative 3" << endl;
        }
        x1 = 0.0; y1 = 0.0; x2 = 2.0; y2 = 1.0; expRes = 0.5;
        der1 = test.FirstDerivative(x1, y1, x2, y2);
        ++nTests;
        if(fabs(der1 - expRes) > 0.0001) {
            ++nFailed;
            cerr << "Error - First derivative 4" << endl;
        }
        x1 = 0.0; y1 = -1.0; x2 = 1.0; y2 = 1.0; expRes = 2.0;
        der1 = test.FirstDerivative(x1, y1, x2, y2);
        ++nTests;
        if(fabs(der1 - expRes) > 0.0001) {
            ++nFailed;
            cerr << "Error - First derivative 5" << endl;
        }
        x1 = 0.0; y1 = 1.0; x2 = 1.0; y2 = 1.0; expRes = 0.0;
        der1 = test.FirstDerivative(x1, y1, x2, y2);
        ++nTests;
        if(fabs(der1 - expRes) > 0.0001) {
            ++nFailed;
            cerr << "Error - First derivative 6" << endl;
        }

        // Second derivative
        x1 = 0.0; y1 = 0.0; x2 = 1.0; y2 = 0.0; x3 = 2.0; y3 = 0.0; expRes = 0.0;
        der2 = test.SecondDerivative(x1, y1, x2, y2, x3, y3);
        ++nTests;
        if(fabs(der2 - expRes) > 0.0001) {
            ++nFailed;
            cerr << "Error - Second derivative 1" << endl;
        }
        x1 = 0.0; y1 = 1.0; x2 = 1.0; y2 = 1.0; x3 = 2.0; y3 = 1.0; expRes = 0.0;
        der2 = test.SecondDerivative(x1, y1, x2, y2, x3, y3);
        ++nTests;
        if(fabs(der2 - expRes) > 0.0001) {
            ++nFailed;
            cerr << "Error - Second derivative 2" << endl;
        }
        x1 = 0.0; y1 = 0.0; x2 = 1.0; y2 = 1.0; x3 = 2.0; y3 = 2.0; expRes = 0.0;
        der2 = test.SecondDerivative(x1, y1, x2, y2, x3, y3);
        ++nTests;
        if(fabs(der2 - expRes) > 0.0001) {
            ++nFailed;
            cerr << "Error - Second derivative 3" << endl;
        }
        x1 = 0.0; y1 = 0.0; x2 = 1.0; y2 = 1.0; x3 = 2.0; y3 = 4.0; expRes = 2.0;
        der2 = test.SecondDerivative(x1, y1, x2, y2, x3, y3);
        ++nTests;
        if(fabs(der2 - expRes) > 0.0001) {
            ++nFailed;
            cerr << "Error - Second derivative 4: Result = " << der2 << "; Expected = " << expRes << endl;
        }

        cout << "Derivatives:  Number of tests = " << nTests << ";  Number of failures = " << nFailed << endl;
    }

    // Test 6 - AutoSet:  y = sin(x)
    // Result: values too concentrated in the points with high absolute second derivative
    {
        cout << "\n\nTest 6: Automatic irregular grid:    y = sin(x)" << endl;
        const string funcName = "y = sin(x)";
        PrecompData<float> itp(funcName);
        const float x0 = 0.0f, x1 = 6.28f;
        itp.SetOversampling(1.5f);
        itp.AutoSet(&TestFuncNonLinSin, x0, x1, nValues);
        cerr << "x0 = " << x0 << "  x1 = " << x1 << "  nValues = " << nValues << endl;  //+T+
        std::vector<float> vx, vy;
        itp.Get(vx, vy);
        cout << "Sizes:  x = " << vx.size() << ";  y = " << vy.size() << endl;
        for(size_t i = 0; i < nValues; ++i) {
            cout << i << ":  " << vx[i] << ", " << vy[i] << endl;
        }
        int n = 100;
        cout << "Compare approximation with real sin(x) function (done on " << n << " points):" << endl;
<<<<<<< .merge_file_R6FDfJ
        float error = 0.0, avgErr = 0.0;
		float minErrX = 1.0e20, minErrY = 1.0e20;
		float maxErrX = 0.0,    maxErrY = 0.0;
=======
        float error = 0.0f, avgErr = 0.0f;
        float minErrX = 1.0e20f, minErrY = 1.0e20f;
        float maxErrX = 0.0f,    maxErrY = 0.0f;
>>>>>>> .merge_file_Kn5us0
        float x = x0, step = (x1 - x0)/n;
        for(int i = 0; i < n; ++i) {
            error = fabs(sin(x) - itp.Interpolate(x));
            cout << i << ": \t" << x << ", \t " << sin(x) << ", \t " << itp.Interpolate(x) << ", \t " << error << endl;
            if(error < minErrY) { minErrX = x; minErrY = error; }
            if(error > maxErrY) { maxErrX = x; maxErrY = error; }
            avgErr += error;
            x += step;
        }
        avgErr /= n;
        cout << "Result:  minimum error = [" << minErrX << ", " << minErrY
             <<      "];  maximum error = [" << maxErrX << ", " << maxErrY
             <<      "];  average error = " << avgErr << endl;
    }


}


} // Utilities


int main()
{
    using namespace Utilities;

    PrecompData_test test;

    return 0;
}
