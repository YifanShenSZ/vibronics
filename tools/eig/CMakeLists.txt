cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(eig)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

enable_language(Fortran)

# Cpp-Library
set(CMAKE_PREFIX_PATH ~/Library/Cpp-Library)
find_package(CL REQUIRED)

# vibron
set(CMAKE_PREFIX_PATH ~/Software/Mine/vibronics/library/vibron)
find_package(vibron REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${vibron_CXX_FLAGS}")

add_executable(eig.exe
    source/steqr.f90
    source/spectrum.cpp
    source/wfn.cpp
    source/main.cpp
)

target_link_libraries(eig.exe
    ${vibron_LIBRARIES}
    ${CL_LIBRARIES}
    iomp5 pthread m dl mkl_intel_lp64 mkl_intel_thread mkl_core
)