# PrecompData

Container of pre-calculated values, returned with interpolation/regression. 

Purpose: improve performance avoiding the realtime computation of complex functions.


#### Currently implemented features

- Precompute regular grid, one dimensional functions, passed as function pointers.
- Zero degree (nearest-neighbor/point sampling/Voronoi) interpolation.


#### TODO

- Copy data on GPU memory.
- First degree (linear) interpolation.
- Precompute regular grid, higher dimensional functions, passed as function pointers.
- Precompute irregular grid, n dimensional functions, passed as function pointers.
- Load data from file.
- Second degree (quadratic) interpolation.


#### Requirements

- C++11
- To copy data on GPU memory:
    - Boost:   http://www.boost.org (from 1.61.0)
    - OpenCL:  http://www.opencl.org


#### Building

