/// PrecompData_test.cpp

/** Test the PrecompData class
 */

#include "PrecompData.h"
#include "PrecompData_test.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using std::cout;
using Utilities::PrecompData;

typedef PrecompData<float, float, 2, 1> pcd21;  // f: (X, Y) --> Z

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

float TestFuncNonLinSin(float x) {      //   y = sin(x)
    return sin(x);
}

pcd21::Y TestFunc21(pcd21::X x) {       //   y = sin(x0) + cos(x1)
    return pcd21::Y { sin(x[0]) + cos(x[1]) };
}


namespace Utilities {


template <typename T = float>
bool TestEqAbs(T value, T expected, T tolerance = 0.01f)
{
    if(fabs(value - expected) <= tolerance)
        return true;    // success

    return false;    // failure
}


template <typename T = float>
bool TestEqRel(T value, T expected, T tolerance = 0.01f)
{
    if(fabs((value - expected)/expected) <= tolerance)
        return true;    // success

    return false;    // failure
}




PrecompData_test::PrecompData_test()
{
	using namespace Utilities;

    cout << std::fixed;

    const int nValues = 20;
    float tol = 0.01f;    // tolerance

    // Test - Conversions RtoI
    {
        cout << "\n\nTest: Conversion real --> index: " << flush;
        const string funcName = "TestFunc";
        PrecompData<> itp(funcName);
        const float x0 = 0.0f, x1 = 6.28f;
        itp.Set(&TestFunc, x0, x1, nValues);
        assert(TestEqAbs(itp.RtoI(x0), size_t(0), size_t(0)) && "Test: Conversion real --> index FAILED on first element.");
        assert(TestEqAbs(itp.RtoI((x1 - x0)/2.0), size_t(nValues/2), size_t(0)) && "Test: Conversion real --> index FAILED on the middle element.");
        assert(TestEqAbs(itp.RtoI(x1), size_t(nValues), size_t(0)) && "Test: Conversion real --> index FAILED on last element.");
        cout << " OK" << endl;
    }

    // Test - Conversions ItoR
    {
        cout << "\n\nTest: Conversion index --> real: " << flush;
        const string funcName = "TestFunc";
        PrecompData<> itp(funcName);
        const float x0 = 0.0f, x1 = 6.28f;
        itp.Set(&TestFunc, x0, x1, nValues);
        assert(TestEqAbs(itp.ItoR(0), x0, tol) && "Test: Conversion index --> real FAILED on first element.");
        assert(TestEqAbs(itp.ItoR(nValues/2), (x1 - x0)/2.0f, tol) && "Test: Conversion index --> real FAILED on the middle element.");
        assert(TestEqAbs(itp.ItoR(nValues), x1, tol) && "Test: Conversion index --> real FAILED on last element.");
        cout << " OK" << endl;
    }

    // Test - Zero-degree interpolation (R --> R)
	{
        cout << "\n\nTest: Zero-degree (nearest-neighbor/point sampling/Voronoi) interpolation:" << endl;
		const string funcName = "TestFunc";
		PrecompData<> itp(funcName);    // default: float type
		itp.SetComment("TestFunc approximation");
		const float x0 = 0.0f, x1 = 6.28f;
		const float step = 0.5f*(x1 - x0)/nValues;
		itp.Set(&TestFunc, x0, x1, nValues);
		float x = x0;
        float err = 0.0f;
        itp.Interpolation(0);
        cout << "Interpolation: " << itp.Interpolation() << endl;
        for(int i = 0; i < nValues; ++i) {
            const float y = itp(x);
            err += fabs(TestFunc(x) - y);
			cout << i << ":\t" << funcName << "(" << x << ") = " << TestFunc(x) << " ~ " << y << endl;
			x += step;
		}
        cout << "Total error = " << err << endl;
	}

    // Test - Linear interpolation (R --> R)
    {
        cout << "\n\nTest: Linear interpolation:" << endl;
        const string funcName = "TestFunc";
        PrecompData<float> itp(funcName);
        const float x0 = 0.0f, x1 = 6.28f;
        const float step = 0.5f*(x1 - x0)/nValues;
        itp.Set(&TestFunc, x0, x1, nValues);
        float x = x0;
        float err = 0.0f;
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

    // Test - AutoSet (R --> R):  y = 2x
    {
        cout << "\n\nTest: Automatic irregular grid:    y = 2x" << endl;
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

    // Test - AutoSet (R --> R):  y = 1/(|x-2| + 0.1)
    {
        cout << "\n\nTest: Automatic irregular grid:    y = 1/(|x-2| + 0.1)" << endl;
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
#if 0
    // Test - Derivatives
    {
        cout << "\n\nTest: Derivatives" << endl;
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

    // Test - AutoSet:  y = sin(x)
    // Result: values too concentrated in the points with high absolute second derivative
    {
        cout << "\n\nTest: Automatic irregular grid:    y = sin(x)" << endl;
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
        float error = 0.0f, avgErr = 0.0f;
        float minErrX = 1.0e20f, minErrY = 1.0e20f;
        float maxErrX = 0.0f,    maxErrY = 0.0f;
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

    // Test - Multidimensions
    {
        cout << "\n\nTest: Storage of data in an NxM space:" << endl;
        const string funcName = "Multidimensions";
        pcd21 itp(funcName);
        itp.SetComment("Y = f(X)    X = x(i,j), Y = y(i)");
        pcd21::X x0 = { 0.00f, 0.00f };     // coordinates of the starting point
        pcd21::X x1 = { 6.28f, 6.28f };     // coordinates of the end point
		const pcd21::X step = { 0.5f*(x1[0] - x0[0])/nValues,
		                        0.5f*(x1[1] - x0[1])/nValues };
        itp.Set(&TestFunc21, x0, x1, nValues*nValues);
        pcd21::X x = x0;
        float err = 0.0f;
        for(int j = 0; j < nValues; ++j)
        {
            x[0] = x0[0];
            for(int i = 0; i < nValues; ++i)
            {
                const pcd21::Y y = itp(x);
                err += fabs(TestFunc21(x)[0] - y[0]);
                cout << i << ":\t" << funcName << "[" << x[0] << ", " << x[1] << "] = " << TestFunc21(x)[0] << " ~ " << y[0] << endl;
                x[0] += step[0];
            }
            x[1] += step[1];
        }
        cout << "Total error = " << err << endl;
    }
#endif
}


} // Utilities


int main()
{
    using namespace Utilities;

    PrecompData_test test;

    return 0;
}
