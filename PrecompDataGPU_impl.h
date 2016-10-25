/**  PrecompDataGPU_impl.h

    Copyright 2016 Pietro Mele
    Apache License 2.0
*/

#ifdef PRECOMPDATA_DEVICE

#include "PrecompData.h"

namespace Utilities {


template<typename T>
int PrecompData<T>::InitDevice()
{
    //+D? device = boost::compute::system::default_device();
    
    return 0;
}


//+TEST
template<typename T>
int PrecompData<T>::CopyOnDevice(boost::compute::device         &device,
                                 boost::compute::context        &context,
                                 boost::compute::command_queue  &queue,
                                 boost::compute::vector<T>      &device_line)
{
    if(line.empty())
        return err_no_data;

    // Resize vector on device
    device_line.resize(line.size(), queue);

    // Copy data from host to device
    boost::compute::copy(
        line.data(), &line.back(), device_line.begin(), queue
    );

    return 0;
}

//+TEST
template<typename T>
int PrecompData<T>::CopyOnDevice(boost::compute::device         **device,
                                 boost::compute::context        **context,
                                 boost::compute::command_queue  **queue,
                                 boost::compute::vector<T>      **device_line)
{
    if(line.empty())
        return err_no_data;
    
    // Create what does not exist, only

    if(*device == nullptr)
    {
        // New device
        *device = compute::system::default_device();

        // New device; previous context and command queue meaningless
        (*context)->~context();
        (*queue)->~command_queue();

        // New context and command queue
        *context = new compute::context(**device);
        *queue   = new compute::command_queue(**context, **device);
    }
    else
    {
        if(*context == nullptr)
        {
            // New context; previous command queue meaningless
            (*queue)->~command_queue();

            // New context and command queue
            *context = new compute::context(**device);
            *queue   = new compute::command_queue(**context, **device);
        }
        else
        {
            if(*queue == nullptr)
            {
                // New command queue
                *queue = new compute::command_queue(**context, **device);
            }
        }
    }
    
    // Resize vector on device
    device_line->resize(line.size(), &queue);

    // Copy data from host to device
    boost::compute::copy(
        line.data(), &line.back(), device_line->begin(), &queue
    );

    return 0;
}


} // Utilities


#else // PRECOMPDATA_DEVICE


#endif // PRECOMPDATA_DEVICE
