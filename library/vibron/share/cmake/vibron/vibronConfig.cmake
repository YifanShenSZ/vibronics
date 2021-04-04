# Find vibron
# -------
#
# Finds vibron
#
# This will define the following variables:
#
#   vibron_FOUND        -- True if the system has vibron
#   vibron_INCLUDE_DIRS -- The include directories for vibron
#   vibron_LIBRARIES    -- Libraries to link against
#
# and the following imported targets:
#
#   vibron

# Find vibron root
# Assume we are in ${vibronROOT}/share/cmake/vibron/vibronConfig.cmake
get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(vibronROOT "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

# include directory
set(vibron_INCLUDE_DIRS ${vibronROOT}/include)

# library
add_library(vibron STATIC IMPORTED)
set(vibron_LIBRARIES vibron)

# dependency 3: Cpp-Library
if(NOT CL_FOUND)
    find_package(CL REQUIRED PATHS ~/Library/Cpp-Library)
    list(APPEND vibron_INCLUDE_DIRS ${CL_INCLUDE_DIRS})
    list(APPEND vibron_LIBRARIES ${CL_LIBRARIES})
endif()

# dependency 2: libtorch
if(NOT TORCH_FOUND)
    find_package(Torch REQUIRED PATHS ~/Software/Programming/libtorch-cuda10.1-1.7.1) 
    list(APPEND vibron_INCLUDE_DIRS ${TORCH_INCLUDE_DIRS})
    list(APPEND vibron_LIBRARIES ${TORCH_LIBRARIES})
    set(vibron_CXX_FLAGS "${vibron_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
endif()

# dependency 1: openmp
if(NOT OPENMP_FOUND)
    find_package(OpenMP REQUIRED)
    set(vibron_CXX_FLAGS "${OpenMP_CXX_FLAGS}")
endif()

# import location
find_library(vibron_LIBRARY vibron PATHS "${vibronROOT}/lib")
set_target_properties(vibron PROPERTIES
    IMPORTED_LOCATION "${vibron_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${vibron_INCLUDE_DIRS}"
    CXX_STANDARD 14
)