/// PrecompData_test.cpp

/** Test the PrecompData class
 */

#include "PrecompData.h"
#include "PrecompData_test.h"
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using Utilities::PrecompData;

typedef PrecompData<100, float, float> pcd;  // f: x --> y

float TestFunc(float x) {
	return sin(x);
}

float TestFuncLin(float x) {            //   y = 2x
    return 2*x;
}

float TestFuncNonLin1(float x) {        //   y = |x|
    return fabs(x);
}

float TestSpike1(float x) {             //   y = 1/(|x-2| + 0.1)   # Spike for x = 2
    return 1/(fabs(x - 2.0f) + 0.1f);
}

float TestFuncNonLinSin(float x) {      //   y = sin(x)
    return sin(x);
}

float TestFuncExpensive(float x) {		// a computationally expensive function
	return (sin(x) + cos(x/2)*sqrt(fabs(x)))/log(fabs(x) + 2.0) + pow(x, x/10.0);
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

	cout << fixed;

	const int    nValues = 20;
	const float  tol     = 0.01f;    // tolerance
	int          verbose = 2;

    // Test - Conversions ScalarToIndex
    {
		cout << "\n\nTest: Conversion scalar --> index: " << flush;
        const string funcName = "TestFunc";
		PrecompData<nValues> itp(funcName);
        const float x0 = 0.0f, x1 = 6.28f;
		itp.set(&TestFunc, x0, x1);
		assert(TestEqAbs(itp.ScalarToIndex(x0), size_t(0), size_t(0)) && "Test: Conversion scalar --> index FAILED on first element.");
		assert(TestEqAbs(itp.ScalarToIndex((x1 - x0)/2.0f), size_t(nValues/2), size_t(0)) && "Test: Conversion scalar --> index FAILED on the middle element.");
		assert(TestEqAbs(itp.ScalarToIndex(x1), size_t(nValues), size_t(0)) && "Test: Conversion scalar --> index FAILED on last element.");
        cout << " OK" << endl;
    }

    // Test - Conversions IndexToScalar
    {
		cout << "\n\nTest: Conversion index --> scalar: " << flush;
        const string funcName = "TestFunc";
		PrecompData<nValues> itp(funcName);
        const float x0 = 0.0f, x1 = 6.28f;
		itp.set(&TestFunc, x0, x1);
		assert(TestEqAbs(itp.IndexToScalar(0), x0, tol) && "Test: Conversion index --> scalar FAILED on first element.");
		assert(TestEqAbs(itp.IndexToScalar(nValues/2), (x1 - x0)/2.0f, tol) && "Test: Conversion index --> scalar FAILED on the middle element.");
		assert(TestEqAbs(itp.IndexToScalar(nValues), x1, tol) && "Test: Conversion index --> scalar FAILED on last element.");
        cout << " OK" << endl;
    }

	// Test - Test mathematical functions
	{
		cout << "\n\nTest: Test mathematical functions: " << endl;
		float x;
		float y, y_ok;

		x = 0.0f;
		y_ok = 0.47619f;
		y = TestSpike1(x);
		cerr << "Expected result = " << y_ok << ";  Actual result = " << y << endl;
		assert(abs(y - y_ok) < 1.0e-2f);

		x = 1.0f;
		y_ok = 0.90909f;
		y = TestSpike1(x);
		cerr << "Expected result = " << y_ok << ";  Actual result = " << y << endl;
		assert(abs(y - y_ok) < 1.0e-2f);

		x = 2.0f;
		y_ok = 10.0f;
		y = TestSpike1(x);
		cerr << "Expected result = " << y_ok << ";  Actual result = " << y << endl;
		assert(abs(y - y_ok) < 1.0e-2f);

		cout << " OK" << endl;
	}

	//+D? Test - Conversions VectorToIndex //+TODO
	/*{
		cout << "\n\nTest: Conversion vector --> index: " << flush;
		const string funcName = "TestFunc";
		pcd12 itp(funcName);
		itp.SetComment("Y = f(X)    X = x(i,j), Y = y(i)");
		const float x0 = 0.00f;
		const float x1 = 6.28f;
		itp.Set(&TestFunc12, x0, x1);

		//itp.Dump();
		itp.Dump(10);
		itp.Dump(-10);

		pcd12::YData y, y_ok;
		y = itp(x0);
		y_ok = TestFunc12(x0);
		cerr << "Expected result = " << y_ok[0] << ";  Actual result = " << y[0] << endl;
		assert(abs(y[0] - y_ok[0]) < 1.0e-2f);

		y = itp(x1);
		y_ok = TestFunc12(x1);
		cerr << "Expected result = " << y_ok[0] << ";  Actual result = " << y[0] << endl;
		assert(abs(y[0] - y_ok[0]) < 1.0e-2f);

		{//+TEMP
			cerr << endl;
			cerr << itp.VectorToIndex(x0) << " = " << size_t(0) << endl;
			cerr << itp.VectorToIndex(x1) << " = " << size_t(nValues*nValues) << endl;
			cerr << itp.VectorToIndex((x1 - x0)/2.0f) << " = " << size_t(nValues*nValues/2) << endl;
		}

		assert(TestEqAbs(itp.VectorToIndex(x0), size_t(0), size_t(0)) && "Test: Conversion vector --> index FAILED on first element.");
		assert(TestEqAbs(itp.VectorToIndex(x1), size_t(nValues*nValues), size_t(0)) && "Test: Conversion vector --> index FAILED on last element.");
		assert(TestEqAbs(itp.VectorToIndex((x1 - x0)/2.0f), size_t(nValues*nValues/2), size_t(0)) && "Test: Conversion vector --> index FAILED on the middle element.");
		cout << " OK" << endl;
	}*/

	// Test - Zero-degree interpolation (R --> R)
	{
        cout << "\n\nTest: Zero-degree (nearest-neighbor/point sampling/Voronoi) interpolation:" << endl;
		const string funcName = "TestFunc";
		PrecompData<nValues> itp(funcName);    // default: float type
		itp.SetComment("TestFunc approximation");
		const float x0 = 0.0f, x1 = 6.28f;
		const float step = 0.5f*(x1 - x0)/nValues;
		itp.set(&TestFunc, x0, x1);
		float x = x0;
        float err = 0.0f;
        itp.Interpolation(0);
        cout << "Interpolation: " << itp.Interpolation() << endl;
        for(int i = 0; i < nValues; ++i) {
            const float y = itp(x);
            err += fabs(TestFunc(x) - y);
			if(verbose > 1)
				cout << i << ":\t" << funcName << "(" << x << ") = " << TestFunc(x) << " ~ " << y << endl;
			x += step;
		}
        cout << "Total error = " << err << endl;
	}

    // Test - Linear interpolation (R --> R)
    {
        cout << "\n\nTest: Linear interpolation:" << endl;
        const string funcName = "TestFunc";
		PrecompData<nValues, float> itp(funcName);
        const float x0 = 0.0f, x1 = 6.28f;
        const float step = 0.5f*(x1 - x0)/nValues;
		itp.set(&TestFunc, x0, x1);
        float x = x0;
        float err = 0.0f;
        itp.Interpolation(1);
        cout << "Interpolation: " << itp.Interpolation() << endl;
        for(int i = 0; i < nValues; ++i) {
            const float y = itp.Interpolate(x);
            err += fabs(TestFunc(x) - y);
			if(verbose > 1)
				cout << i << ":\t" << funcName << "(" << x << ") = " << TestFunc(x) << " ~ " << y << endl;
            x += step;
        }
		cout << "Total error = " << err << endl;
    }

    // Test - AutoSet (R --> R):  y = 2x
    {
        const string funcName = "y = 2x";
		cout << "\n\nTest: Automatic irregular grid:    " << funcName << endl;
		PrecompData<nValues, float> itp(funcName);
        const float x0 = 0.0f, x1 = 6.28f;
		itp.AutoSet(&TestFuncLin, x0, x1);
        cerr << "x0 = " << x0 << "  x1 = " << x1 << "  nValues = " << nValues << endl;  //+T+
        std::vector<float> vx, vy;
		itp.get(vx, vy);
        cout << "Sizes:  x = " << vx.size() << ";  y = " << vy.size() << endl;
		if(verbose > 1)
			for(size_t i = 0; i < nValues; ++i) {
				cout << i << ":  " << vx[i] << ", " << vy[i] << endl;
			}
    }

	// Test - Regular grid (R --> R):  y = 1/(|x-2| + 0.1)   # Spike for x = 2
	{
		const string funcName = "y = 1/(|x-2| + 0.1)";
		cout << "\n\nTest: Regular grid:    " << funcName << endl;
		PrecompData<nValues, float> itp(funcName);
		const float x0 = 0.0f, x1 = 6.28f;
		itp.set(&TestSpike1, x0, x1);
		cerr << "x0 = " << x0 << "  x1 = " << x1 << "  nValues = " << nValues << endl;  //+T+
		std::vector<float> vx, vy;
		itp.get(vx, vy);
		cout << "Sizes:  x = " << vx.size() << ";  y = " << vy.size() << endl;
		if(verbose > 1)
			for(size_t i = 0; i < nValues; ++i) {
				cout << i << ":  " << vx[i] << ", " << vy[i] << endl;
			}
		float absErrKnown = itp.EvaluateAbsErrorKnownData();
		float absErr      = itp.EvaluateAbsError(2*nValues);
		cout << "Absolute error:  on known data = " << absErrKnown
		     <<                "; on random data = " << absErr << endl;
	}

	// Test - AutoSet (R --> R):  y = 1/(|x-2| + 0.1)   # Spike for x = 2
	{
		const string funcName = "y = 1/(|x-2| + 0.1)";
		cout << "\n\nTest: Automatic irregular grid:    " << funcName << endl;
		PrecompData<nValues, float> itp(funcName);
		const float x0 = 0.0f, x1 = 6.28f;
		itp.AutoSet(&TestSpike1, x0, x1);
		cerr << "x0 = " << x0 << "  x1 = " << x1 << "  nValues = " << nValues << endl;  //+T+
		std::vector<float> vx, vy;
		itp.get(vx, vy);
		cout << "Sizes:  x = " << vx.size() << ";  y = " << vy.size() << endl;
		if(verbose > 1)
			for(size_t i = 0; i < nValues; ++i) {
				cout << i << ":  " << vx[i] << ", " << vy[i] << endl;
			}
		float absErrKnown = itp.EvaluateAbsErrorKnownData();
		float absErr      = itp.EvaluateAbsError(2*nValues);
		cout << "Absolute error:  on known data = " << absErrKnown
		     <<                "; on random data = " << absErr << endl;
	}

	// Test - Derivatives
    {
        cout << "\n\nTest: Derivatives" << endl;
        int nTests = 0, nFailed = 0;
        const string funcName = "Derivatives";
        float x1, y1, x2, y2, x3, y3, der1, der2, expRes;
		PrecompData<nValues, float> test;

        // First derivative
        x1 = 0.0; y1 = 0.0; x2 = 1.0; y2 = 0.0; expRes = 0.0;
        der1 = test.FirstDerivative(x1, y1, x2, y2);
        ++nTests;
		if(fabs(der1 - expRes) > 0.0001f) {
            ++nFailed;
            cerr << "Error - First derivative 1" << endl;
        }
        x1 = 0.0; y1 = 0.0; x2 = 1.0; y2 = 1.0; expRes = 1.0;
        der1 = test.FirstDerivative(x1, y1, x2, y2);
        ++nTests;
		if(fabs(der1 - expRes) > 0.0001f) {
            ++nFailed;
            cerr << "Error - First derivative 2" << endl;
        }
        x1 = 1.0; y1 = 0.0; x2 = 0.0; y2 = 1.0; expRes = -1.0;
        der1 = test.FirstDerivative(x1, y1, x2, y2);
        ++nTests;
		if(fabs(der1 - expRes) > 0.0001f) {
            ++nFailed;
            cerr << "Error - First derivative 3" << endl;
        }
        x1 = 0.0; y1 = 0.0; x2 = 2.0; y2 = 1.0; expRes = 0.5;
        der1 = test.FirstDerivative(x1, y1, x2, y2);
        ++nTests;
		if(fabs(der1 - expRes) > 0.0001f) {
            ++nFailed;
            cerr << "Error - First derivative 4" << endl;
        }
        x1 = 0.0; y1 = -1.0; x2 = 1.0; y2 = 1.0; expRes = 2.0;
        der1 = test.FirstDerivative(x1, y1, x2, y2);
        ++nTests;
		if(fabs(der1 - expRes) > 0.0001f) {
            ++nFailed;
            cerr << "Error - First derivative 5" << endl;
        }
        x1 = 0.0; y1 = 1.0; x2 = 1.0; y2 = 1.0; expRes = 0.0;
        der1 = test.FirstDerivative(x1, y1, x2, y2);
        ++nTests;
		if(fabs(der1 - expRes) > 0.0001f) {
            ++nFailed;
            cerr << "Error - First derivative 6" << endl;
        }

        // Second derivative
        x1 = 0.0; y1 = 0.0; x2 = 1.0; y2 = 0.0; x3 = 2.0; y3 = 0.0; expRes = 0.0;
        der2 = test.SecondDerivative(x1, y1, x2, y2, x3, y3);
        ++nTests;
		if(fabs(der2 - expRes) > 0.0001f) {
            ++nFailed;
            cerr << "Error - Second derivative 1" << endl;
        }
        x1 = 0.0; y1 = 1.0; x2 = 1.0; y2 = 1.0; x3 = 2.0; y3 = 1.0; expRes = 0.0;
        der2 = test.SecondDerivative(x1, y1, x2, y2, x3, y3);
        ++nTests;
		if(fabs(der2 - expRes) > 0.0001f) {
            ++nFailed;
            cerr << "Error - Second derivative 2" << endl;
        }
        x1 = 0.0; y1 = 0.0; x2 = 1.0; y2 = 1.0; x3 = 2.0; y3 = 2.0; expRes = 0.0;
        der2 = test.SecondDerivative(x1, y1, x2, y2, x3, y3);
        ++nTests;
		if(fabs(der2 - expRes) > 0.0001f) {
            ++nFailed;
            cerr << "Error - Second derivative 3" << endl;
        }
        x1 = 0.0; y1 = 0.0; x2 = 1.0; y2 = 1.0; x3 = 2.0; y3 = 4.0; expRes = 2.0;
        der2 = test.SecondDerivative(x1, y1, x2, y2, x3, y3);
        ++nTests;
		if(fabs(der2 - expRes) > 0.0001f) {
            ++nFailed;
            cerr << "Error - Second derivative 4: Result = " << der2 << "; Expected = " << expRes << endl;
        }

        cout << "Derivatives:  Number of tests = " << nTests << ";  Number of failures = " << nFailed << endl;
    }

    // Test - AutoSet:  y = sin(x)
	// Result: values too concentrated in the points with high absolute second derivative.
	//         Large errors at the end of the interval.
    {
        const string funcName = "y = sin(x)";
		cout << "\n\nTest: Automatic irregular grid:    " << funcName << endl;
		PrecompData<nValues, float> itp(funcName);
        const float x0 = 0.0f, x1 = 6.28f;
        itp.SetOversampling(1.5f);
		itp.AutoSet(&TestFuncNonLinSin, x0, x1);
        cerr << "x0 = " << x0 << "  x1 = " << x1 << "  nValues = " << nValues << endl;  //+T+
        std::vector<float> vx, vy;
		itp.get(vx, vy);
        cout << "Sizes:  x = " << vx.size() << ";  y = " << vy.size() << endl;
		if(verbose > 1)
			for(size_t i = 0; i < nValues; ++i) {
				cout << i << ":  " << vx[i] << ", " << vy[i] << endl;
			}
        int n = 100;
		cout << "Compare approximation with real " << funcName << " function (done on " << n << " points):" << endl;
        float error = 0.0f, avgErr = 0.0f;
        float minErrX = 1.0e20f, minErrY = 1.0e20f;
        float maxErrX = 0.0f,    maxErrY = 0.0f;
        float x = x0, step = (x1 - x0)/n;
        for(int i = 0; i < n; ++i) {
			float interp;
			const size_t r = itp.Interpolate(x, interp);
			if(r == PrecompData_base::wrn_x_less_than_min || r == PrecompData_base::wrn_x_more_than_max) { cout << "Out of range" << endl; }
			error = fabs(sin(x) - interp);
			if(verbose > 1)
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

	// Test - Performance of computed vs precomputed function
	{
		cout << "\n\nTest: Performance of computed vs precomputed function:" << endl;
		const string funcName = "TestFunc";
		const int multiplier = 1000;
		PrecompData<nValues, float> itp(funcName);
		const float x0 = 0.0f, x1 = 6.28f;
		const float step = 0.5f*(x1 - x0)/(multiplier*nValues);
		itp.set(&TestFuncExpensive, x0, x1);
		float x = x0;
		float err = 0.0f;
		itp.Interpolation(1);
		cout << "Interpolation: " << itp.Interpolation() << endl;
		chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();
		chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();
		chrono::duration<float> timeComp = start - start;
		chrono::duration<float> timePrecomp = start - start;
		x = x0;
		start = chrono::system_clock::now();
		for(int i = 0; i < multiplier*nValues; ++i)
		{
			const float y_comp = TestFuncExpensive(x);
			x += step;
		}
		end = chrono::system_clock::now();
		timeComp += end - start;
		x = x0;
		start = chrono::system_clock::now();
		for(int i = 0; i < multiplier*nValues; ++i)
		{
			const float y_prec = itp.Interpolate(x);
			x += step;
		}
		end = chrono::system_clock::now();
		timePrecomp += end - start;
		x = x0;
		for(int i = 0; i < multiplier*nValues; ++i)
		{
			const float y_comp = TestFuncExpensive(x);
			const float y_prec = itp.Interpolate(x);
			err += fabs(y_prec - y_comp);
			x += step;
		}
		typedef chrono::duration<float> float_seconds;
		auto timeCompSecs = std::chrono::duration_cast<float_seconds>(timeComp);
		auto timePrecompSecs = std::chrono::duration_cast<float_seconds>(timePrecomp);
		cout << "Total error = " << err;
		cout << "\nCompute time = " << timeCompSecs.count();
		cout << "\nPreComp time = " << timePrecompSecs.count();
		cout << "\nImprovement  = " << 100.0*(timeComp - timePrecomp)/timeComp << "%" << endl;
	}

	// Test - Performance of computed vs precomputed function (linear interpolation)
	{
		cout << "\n\nTest: Performance of computed vs precomputed function:" << endl;
		const string funcName = "TestFunc";
		PrecompData<nValues, float> itp(funcName);
		const float x0 = 0.0f, x1 = 6.28f;
		const int nTestPoints = 10000;
		itp.set(&TestFuncExpensive, x0, x1);
		itp.Interpolation(1);
		cout << "Interpolation: " << itp.Interpolation() << endl;
		cout << "Performance improvement = " << itp.PerformanceImprovement(nTestPoints) << endl;
		cout << "Total error             = " << itp.EvaluateAbsError(10*nValues) << endl;
	}

	cout << endl;
}


} // Utilities


int main()
{
    using namespace Utilities;

    PrecompData_test test;

    return 0;
}
