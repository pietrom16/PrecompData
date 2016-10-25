# PrecompData - TODO

#### Copy data on device memory.

- How to make 'context' and 'queue' member variables?
- How to share 'device_line' with the other tools? Getting 'device_line' as a parameter?

- Anything like:

    int CopyOnDevice(const compute::device&, 
                     const compute::context&, 
                     const compute::command_queue&, 
                     const compute::vector<T> &device_line);
