# PrecompData

Container of pre-calculated values, returned with interpolation/regression.

Purpose: improve performance avoiding the realtime computation of complex functions.


### Currently implemented features

- Get data through function pointers.
- Precompute regular grid, one dimensional functions.
- Zero degree (nearest-neighbor/point sampling/Voronoi) interpolation.
- First degree (linear) interpolation for functions of one variable.
- Copy data on GPU/device memory. [TEST]
- Precompute irregular grid, one dimensional functions. [OPTIMIZE/TEST]
- Test whether it is worth to precompute the data, in terms of performance, approximation error and memory requirements.


### TODO

- Precompute regular grid, n dimensional functions.
- Precompute irregular grid, n dimensional functions.
- Load data from file.
- Second degree (quadratic) interpolation.
- Spline interpolation.


### Requirements

- C++11
- To copy data on GPU/device memory:
	- Boost:   http://www.boost.org (from 1.61.0)
	- OpenCL:  https://www.khronos.org/opencl/


### Building

- To enable coping data on GPU/device memory:
	- `#define PRECOMPDATA_DEVICE` either in the `PrecompData.h` header or in the build system.
	- Add the Boost compute include directory to the compilation flags.
	- Link with the system's OpenCL library.
	- GCC example:  `g++ -I/path/to/compute/include main.cpp -lOpenCL`

### Examples

- Here is a basic example of PrecompData usage:

```C++
#include "PrecompData.h"
#include <cmath>
using namespace std;

// Function to precalculate
float MyFunction(float x) {
    return  sin(x) + cos(2*x);
}

int main()
{
    using namespace Utilities;

    const float x0 = 0.0f, x1 = 6.28f;
    const int nValues = 10;

    PrecompData<float> func;

    func.Set(&MyFunction, x0, x1, nValues);

    float x = 1.234;
    float y = func.Interpolate(x);
}
```

- For a full set of examples, look at the test file `PrecompData_test.cpp`.

### Error reporting

    - The Status() member function returns the latest error/warning/success value:
	    Success: 0
		Error:   [-1, -100]
		Warning: < -100
	- When any member function is called, the status is reset.


### References

- ["Sequential and Parallel Generation of Unstructured 3D Meshes"](http://mech.fsv.cvut.cz/~dr/papers/Thesis98/thesis.html)
    - Chapter: ["Curve Discretization"](http://mech.fsv.cvut.cz/~dr/papers/Thesis98/node40.html)  
- ["Space-filling curve"](https://en.wikipedia.org/wiki/Space-filling_curve)  
