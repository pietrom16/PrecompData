/**  PrecompDataGPU_impl.h

    Copyright 2016 Pietro Mele
    Apache License 2.0
*/

#ifdef PRECOMPDATA_GPU

#include "PrecompData.h"

namespace Utilities {


template<typename T>
int PrecompData<T>::InitGPU()
{
    //+D? device = boost::compute::system::default_device();
    
    return 0;
}


template<typename T>
int PrecompData<T>::CopyOnGPU(boost::compute::device         &device,
                              boost::compute::context        &context,
                              boost::compute::command_queue  &queue,
                              boost::compute::vector<T>      &device_line)
{
    if(line.empty())
        return err_no_data;

    // Resize vector on device
    device_line.resize(line.size(), context);

    // Copy data from host to device
    boost::compute::copy(
        line.data(), &line.back(), device_line.begin(), queue
    );

    return 0;
}


template<typename T>
int PrecompData<T>::CopyOnGPU(boost::compute::device         *device,
                              boost::compute::context        *context,
                              boost::compute::command_queue  *queue,
                              boost::compute::vector<T>      *device_line)
int PrecompData<T>::CopyOnGPU(boost::compute::device         **device,
                              boost::compute::context        **context,
                              boost::compute::command_queue  **queue,
                              boost::compute::vector<T>      **device_line)
{
    //+TODO - Create if null

    if(line.empty())
        return err_no_data;

    // Resize vector on device
    device_line->resize(line.size(), &context);

    // Copy data from host to device
    boost::compute::copy(
        line.data(), &line.back(), device_line->begin(), &queue
    );

    return 0;
}


} // Utilities


#else // PRECOMPDATA_GPU


#endif // PRECOMPDATA_GPU
