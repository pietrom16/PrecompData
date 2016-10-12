/// PrecompData_test.cpp

/** Test the PrecompData class
 */

#include "PrecompData.h"
#include <cmath>
#include <vector>

using namespace std;

float TestFunc(float y) {
	return sin(y);
}


int main()
{
	using namespace Utilities;

	const std::string funcName = "TestFunc";

	PrecompData<float> itp1("TestFunc");
	PrecompData<float> itp2(funcName);
	//itp.SetComment("TestFunc approximation");

	//itp.Set(&TestFunc, 0.0, 6.28, 10);

	return 0;
}

