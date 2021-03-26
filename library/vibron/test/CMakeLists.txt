cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(test)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Release)

# vibron
set(CMAKE_PREFIX_PATH ~/Software/Mine/vibronics/library/vibron)
find_package(vibron REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${vibron_CXX_FLAGS}")

add_executable(test.exe main.cpp)

target_link_libraries(test.exe "${vibron_LIBRARIES}")
