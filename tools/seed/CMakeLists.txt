cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Release)

# Cpp-Library
set(CMAKE_PREFIX_PATH ~/Library/Cpp-Library)
find_package(CL REQUIRED)

# vibron
set(CMAKE_PREFIX_PATH ~/Software/Mine/vibronics/library/vibron)
find_package(vibron REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${vibron_CXX_FLAGS}")

# seed
set(CMAKE_PREFIX_PATH ~/Software/Mine/vibronics/library/seed)
find_package(seed REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${seed_CXX_FLAGS}")

add_executable(seed.exe main.cpp)

target_link_libraries(seed.exe
    ${seed_LIBRARIES}
    ${vibron_LIBRARIES}
    ${CL_LIBRARIES}
)