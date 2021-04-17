cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(v2)
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

# Lanczos
set(CMAKE_PREFIX_PATH ~/Software/Mine/vibronics/library/Lanczos)
find_package(Lanczos REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Lanczos_CXX_FLAGS}")

# seed
set(CMAKE_PREFIX_PATH ~/Software/Mine/vibronics/library/seed)
find_package(seed REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${seed_CXX_FLAGS}")

add_executable(v2.exe main.cpp)

target_link_libraries(v2.exe
    ${seed_LIBRARIES} ${Lanczos_LIBRARIES}
    ${vibron_LIBRARIES} ${CL_LIBRARIES}
)