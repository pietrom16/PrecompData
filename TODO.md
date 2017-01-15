# PrecompData - TODO


- Making transition to multidimensional data...
    - Precompute regular grid, n dimensional functions.
	- Precompute irregular grid, one dimensional functions. [OPTIMIZE/TEST]
	- Precompute irregular grid, n dimensional functions.
- Load data from file.
- Second degree (quadratic) interpolation.
- Test whether it is worth to precompute the data, in terms of performance and memory requirements.
- Add a version member variable.
- Check input parameters of most functions (assert).
- Consider using a proper matrix library.
- Option to add noise to data.
- Dynamic resolution, i.e. add new points to the pre-calculated ones when they are in regions not currently covered.


#### Optimizations

- Check for eventual performance penalties using boost::multi_array vs std::vector vs std::array for one-dimensional data.


##### Copy data on device memory.

- Fix device related linker error: Error	LNK1104	cannot open file 'C:\libs\boost\boost_1_62_0\stage\lib\.obj'
- How to make 'context' and 'queue' member variables?
- How to share 'device_line' with the other tools? Getting 'device_line' as a parameter?


#### References

["SO/templates"](http://stackoverflow.com/questions/610245/where-and-why-do-i-have-to-put-the-template-and-typename-keywords)
http://www.boost.org/doc/libs/1_62_0/libs/multi_array/doc/index.html
https://arxiv.org/pdf/1608.00099.pdf
https://bitbucket.org/orserang/triot/
