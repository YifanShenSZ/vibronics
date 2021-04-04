# Find seed
# -------
#
# Finds seed
#
# This will define the following variables:
#
#   seed_FOUND        -- True if the system has seed
#   seed_INCLUDE_DIRS -- The include directories for seed
#   seed_LIBRARIES    -- Libraries to link against
#
# and the following imported targets:
#
#   seed

# Find seed root
# Assume we are in ${seedROOT}/share/cmake/seed/seedConfig.cmake
get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(seedROOT "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

# include directory
set(seed_INCLUDE_DIRS ${seedROOT}/include)

# library
add_library(seed STATIC IMPORTED)
set(seed_LIBRARIES seed)

# dependency 3: vibron
if(NOT vibron_FOUND)
    find_package(vibron REQUIRED PATHS ~/Software/Mine/vibronics/library/vibron)
    list(APPEND seed_INCLUDE_DIRS ${vibron_INCLUDE_DIRS})
    list(APPEND seed_LIBRARIES ${vibron_LIBRARIES})
    set(seed_CXX_FLAGS "${seed_CXX_FLAGS} ${vibron_CXX_FLAGS}")
endif()

# dependency 2: Cpp-Library
if(NOT CL_FOUND)
    find_package(CL REQUIRED PATHS ~/Library/Cpp-Library)
    list(APPEND seed_INCLUDE_DIRS ${CL_INCLUDE_DIRS})
    list(APPEND seed_LIBRARIES ${CL_LIBRARIES})
endif()

# dependency 1: openmp
if(NOT OPENMP_FOUND)
    find_package(OpenMP REQUIRED)
    set(seed_CXX_FLAGS "${OpenMP_CXX_FLAGS}")
endif()

# import location
find_library(seed_LIBRARY seed PATHS "${seedROOT}/lib")
set_target_properties(seed PROPERTIES
    IMPORTED_LOCATION "${seed_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${seed_INCLUDE_DIRS}"
    CXX_STANDARD 14
)