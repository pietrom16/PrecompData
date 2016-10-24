# PrecompData

Container of pre-calculated values, returned with interpolation/regression. 

Purpose: improve performance avoiding the realtime computation of complex functions.


#### Currently implemented features

- Precompute regular grid, one dimensional functions, passed as function pointers.
- Zero degree interpolation.


#### TODO

- First degree (linear) interpolation.
- Copy data on GPU memory.
- Precompute regular grid, higher dimensional functions, passed as function pointers.
- Precompute irregular grid, n dimensional functions, passed as function pointers.
- Load data from file.
- Second degree (quadratic) interpolation.

