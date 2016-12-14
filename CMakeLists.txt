# CMakeLists.txt - PrecompData
#
# Copyright 2016 Pietro Mele

cmake_minimum_required (VERSION 3.0)

set(PRJ PrecompData)

message("==============================================")
message("Building project: ${PRJ}")
message("==============================================")

project(${PRJ})

# Version number
set (${PRJ}_VERSION_MAJOR 1)
set (${PRJ}_VERSION_MINOR 3)
set (${PRJ}_VERSION_PATCH 0)

# add_definitions (-DPRECOMPDATA_DEVICE)

# Specify the build type on the command line: Release, Debug, RelWithDebInfo, MinSizeRel
if(CMAKE_BUILD_TYPE MATCHES "")
	set(CMAKE_BUILD_TYPE Debug)
	#set(CMAKE_BUILD_TYPE Release)
endif()

# Specify C++ standard version
set(CMAKE_CXX_FLAGS "-std=c++14 -stdlib=libc++")

# Extra code in debug mode
if(CMAKE_BUILD_TYPE MATCHES "Debug")
	message("Building extra checks and output in debug mode.")
	add_definitions(-DVERBOSE)
endif()


message("Building with: " ${CMAKE_CXX_COMPILER} " " ${CMAKE_CXX_FLAGS} " " ${CMAKE_BUILD_TYPE})

set(SRC
	PrecompData_base.h
	PrecompData_base.cpp
	PrecompData.h
	PrecompData_impl.h
	PrecompData_multi.h
	PrecompData_multi_impl.h
	PrecompDataDevice_impl.h
	PrecompData_dump.h
	PrecompData_test.h
	PrecompData_test.cpp
	StlExt.hpp
	README.md
	TODO.md
)

add_executable(${PRJ}_test ${SRC})


# Out of source build:
#    cd <build_dir>
#    cmake <path_to_source>
