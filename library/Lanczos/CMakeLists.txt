cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Release)

# OpenMP
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Cpp-Library
set(CMAKE_PREFIX_PATH ~/Library/Cpp-Library)
find_package(CL REQUIRED)

# vibron
set(CMAKE_PREFIX_PATH ~/Software/Mine/vibronics/library/vibron)
find_package(vibron REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${vibron_CXX_FLAGS}")

include_directories(include)

add_library(Lanczos STATIC
    source/Hamiltonian.cpp
    source/mv.cpp
    source/Lanczos.cpp
)

target_link_libraries(Lanczos "${vibron_LIBRARIES}" "${CL_LIBRARIES}")