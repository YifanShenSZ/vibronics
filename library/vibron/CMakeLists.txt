cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Release)

# OpenMP
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# libtorch
set(CMAKE_PREFIX_PATH ~/Software/Programming/libtorch-cuda10.1-1.7.1)
find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

# Cpp-Library
set(CMAKE_PREFIX_PATH ~/Library/Cpp-Library)
find_package(CL REQUIRED)

include_directories(include)

add_library(vibron STATIC
    source/harmonic.cpp
    source/vibration/vibration.cpp
    source/vibration/set.cpp
    source/options.cpp
    source/wfn.cpp
)

target_link_libraries(vibron ${CL_LIBRARIES} ${TORCH_LIBRARIES})