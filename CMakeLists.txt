cmake_minimum_required (VERSION 3.0)
project (PrecompData)

set (PROJECT PrecompData)

# Version number
set (${PROJECT}_VERSION_MAJOR 1)
set (${PROJECT}_VERSION_MINOR 1)
set (${PROJECT}_VERSION_PATCH 0)


add_executable(PrecompData_test PrecompData_test.cpp)

# Out of source build:
#    cd <build_dir>
#    cmake <path_to_source>
