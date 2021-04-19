# Find plot
# -------
#
# Finds plot
#
# This will define the following variables:
#
#   plot_FOUND        -- True if the system has plot
#   plot_INCLUDE_DIRS -- The include directories for plot
#   plot_LIBRARIES    -- Libraries to link against
#
# and the following imported targets:
#
#   plot

# Find plot root
# Assume we are in ${plotROOT}/share/cmake/plot/plotConfig.cmake
get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(plotROOT "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

# include directory
set(plot_INCLUDE_DIRS ${plotROOT}/include)

# library
add_library(plot STATIC IMPORTED)
set(plot_LIBRARIES plot)

# dependency 3: vibron
if(NOT vibron_FOUND)
    find_package(vibron REQUIRED PATHS ~/Software/Mine/vibronics/library/vibron)
    list(APPEND plot_INCLUDE_DIRS ${vibron_INCLUDE_DIRS})
    list(APPEND plot_LIBRARIES ${vibron_LIBRARIES})
    set(plot_CXX_FLAGS "${plot_CXX_FLAGS} ${vibron_CXX_FLAGS}")
endif()

# dependency 2: Cpp-Library
if(NOT CL_FOUND)
    find_package(CL REQUIRED PATHS ~/Library/Cpp-Library)
    list(APPEND plot_INCLUDE_DIRS ${CL_INCLUDE_DIRS})
    list(APPEND plot_LIBRARIES ${CL_LIBRARIES})
endif()

# dependency 1: openmp
if(NOT OPENMP_FOUND)
    find_package(OpenMP REQUIRED)
    set(plot_CXX_FLAGS "${OpenMP_CXX_FLAGS}")
endif()

# import location
find_library(plot_LIBRARY plot PATHS "${plotROOT}/lib")
set_target_properties(plot PROPERTIES
    IMPORTED_LOCATION "${plot_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${plot_INCLUDE_DIRS}"
    CXX_STANDARD 14
)