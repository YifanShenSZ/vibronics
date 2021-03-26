# Find Lanczos
# -------
#
# Finds Lanczos
#
# This will define the following variables:
#
#   Lanczos_FOUND        -- True if the system has Lanczos
#   Lanczos_INCLUDE_DIRS -- The include directories for Lanczos
#   Lanczos_LIBRARIES    -- Libraries to link against
#
# and the following imported targets:
#
#   Lanczos

# Find Lanczos root
# Assume we are in ${LanczosROOT}/share/cmake/Lanczos/LanczosConfig.cmake
get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(LanczosROOT "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

# include directory
set(Lanczos_INCLUDE_DIRS ${LanczosROOT}/include)

# library
add_library(Lanczos STATIC IMPORTED)
set(Lanczos_LIBRARIES Lanczos)

# dependency 1: openmp
if(NOT OPENMP_FOUND)
    find_package(OpenMP REQUIRED)
    set(Lanczos_CXX_FLAGS "${OpenMP_CXX_FLAGS}")
endif()

# dependency 2: Cpp-Library
if(NOT CL_FOUND)
    find_package(CL REQUIRED PATHS ~/Library/Cpp-Library)
    list(APPEND Lanczos_INCLUDE_DIRS ${CL_INCLUDE_DIRS})
    list(APPEND Lanczos_LIBRARIES ${CL_LIBRARIES})
endif()

# dependency 3: vibron
if(NOT vibron_FOUND)
    find_package(vibron REQUIRED PATHS ~/Software/Mine/vibronics/library/vibron)
    list(APPEND Lanczos_INCLUDE_DIRS ${vibron_INCLUDE_DIRS})
    list(APPEND Lanczos_LIBRARIES ${vibron_LIBRARIES})
    set(Lanczos_CXX_FLAGS "${Lanczos_CXX_FLAGS} ${vibron_CXX_FLAGS}")
endif()

# import location
find_library(Lanczos_LIBRARY Lanczos PATHS "${LanczosROOT}/lib")
set_target_properties(Lanczos PROPERTIES
    IMPORTED_LOCATION "${Lanczos_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${Lanczos_INCLUDE_DIRS}"
    CXX_STANDARD 14
)