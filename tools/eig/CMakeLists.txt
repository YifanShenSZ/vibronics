cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(eig)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

enable_language(Fortran)

# Cpp-Library
set(CMAKE_PREFIX_PATH ~/Library/Cpp-Library)
find_package(CL REQUIRED)

add_executable(eig.exe
    source/steqr.f90
    source/main.cpp
)

target_link_libraries(eig.exe
    ${CL_LIBRARIES}
    iomp5 pthread m dl mkl_intel_lp64 mkl_intel_thread mkl_core
)