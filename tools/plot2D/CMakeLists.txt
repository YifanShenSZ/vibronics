cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Release)

# Cpp-Library
set(CMAKE_PREFIX_PATH ~/Library/Cpp-Library)
find_package(CL REQUIRED)

# Torch-Chemistry
set(CMAKE_PREFIX_PATH ~/Library/Torch-Chemistry)
find_package(tchem REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TCHEM_CXX_FLAGS}")

# Lanczos
set(CMAKE_PREFIX_PATH ~/Software/Mine/vibronics/library/Lanczos)
find_package(Lanczos REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${vibron_CXX_FLAGS}")

# plot
set(CMAKE_PREFIX_PATH ~/Software/Mine/vibronics/library/plot)
find_package(plot REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${plot_CXX_FLAGS}")

add_executable(plot2D.exe
    source/HdKernel.cpp
    source/NCKernel.cpp
    source/main.cpp
)

include_directories(include)

target_link_libraries(plot2D.exe
    ${plot_LIBRARIES} ${Lanczos_LIBRARIES}
    ${tchem_LIBRARIES} ${CL_LIBRARIES}
)