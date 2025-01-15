cmake_minimum_required(VERSION 3.14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

project(cfm LANGUAGES CXX)

#It is recommended, but not strictly necessary, to build and run the tests
set(INCLUDE_TESTS ON CACHE BOOL "Do you want to build the tests?")

#Building the training code requires extra dependencies, e.g. MPI
set(INCLUDE_TRAIN ON CACHE BOOL "Do you want to build the training code?")

# BOOST
set(Boost_USE_STATIC_LIBS OFF CACHE BOOL "Use static libraries from Boost")
include(FindBoost)

set(Boost_components filesystem system serialization program_options thread)
# For some reason msvc will give LNK2019 if not include thread

if (INCLUDE_TESTS)
    list(APPEND Boost_components unit_test_framework)
endif()

find_package(Boost REQUIRED COMPONENTS ${Boost_components})
include_directories(${Boost_INCLUDE_DIRS})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


#include cmake finds
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
# RDKit
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
include(FindRDKit)

if (RDKIT_FOUND)
    include_directories(${RDKIT_INCLUDE_DIR} ${RDKIT_INCLUDE_EXT_DIR})
    link_directories(${RDKIT_LIBRARY_DIR})
else ()
    message(STATUS "RDKIT_INCLUDE_DIR ${RDKIT_INCLUDE_DIR}")
    message(STATUS "RDKIT_INCLUDE_DIR ${RDKIT_INCLUDE_EXT_DIR}")
    message(STATUS "RDKIT_INCLUDE_DIR ${RDKIT_LIBRARIES}")
    message(FATAL_ERROR "Could not find RDKit")
endif ()

# LPSolve
set(LPSOLVE_INCLUDE_DIR "" CACHE STRING "Where are the LPSolve headers?")
set(LPSOLVE_LIBRARY_DIR "" CACHE STRING "Where are the LPSolve libraries?")
message(STATUS "Search For lPsolve lib at ${LPSOLVE_LIBRARY_DIR} and headers at ${LPSOLVE_INCLUDE_DIR} ")


find_path(LPSOLVE_LIBRARY_HEADERS lp_lib.h PATHS
                        ${LPSOLVE_INCLUDE_DIR}
                        /usr/include/lpsolve
                        /usr/local/include/lpsolve
                        /opt/lpsolve
                        ~/lpsolve
                        /usr/include/lp_solve
                        /usr/local/include/lp_solve
                        /opt/lp_solve
                        ~/lp_solve)

if (NOT LPSOLVE_LIBRARY_HEADERS)
    MESSAGE(FATAL_ERROR "Cannot find LPSOLVE Headers!")
else()
    MESSAGE(STATUS "Found LPSOLVE Headers at ${LPSOLVE_LIBRARY_HEADERS}")
endif (NOT LPSOLVE_LIBRARY_HEADERS)

include(FetchContent)
FetchContent_Declare(coretools
        GIT_REPOSITORY https://bitbucket.org/wegmannlab/coretools.git
        GIT_TAG madleina/interface_stattools
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/coretools"
)
FetchContent_MakeAvailable(coretools)
target_link_libraries(coretools)

find_library(LPSOLVE_LIBRARY lpsolve55 PATHS
                         ${LPSOLVE_LIBRARY_DIR}
                        /usr/lib
                        /usr/lib/lpsolve
                        /usr/local/lib/
                        /usr/local/lib/lpsolve
                        /opt/lpsolve
                        ~/lpsolve
                        /usr/lib/lp_solve
                        /usr/local/lib/lp_solve
                        /opt/lp_solve
                        ~/lp_solve)

if (NOT LPSOLVE_LIBRARY)
    MESSAGE(FATAL_ERROR "Cannot find LPSOLVE libraries!")
else()
    MESSAGE(STATUS "Found LPSOLVE libraries file ${LPSOLVE_LIBRARY}")
endif (NOT LPSOLVE_LIBRARY)

include_directories(${LPSOLVE_LIBRARY_HEADERS})
link_directories(${LPSOLVE_LIBRARY_HEADERS})
mark_as_advanced(${LPSOLVE_LIBRARY})

find_package(Threads)
find_package(OpenMP)

set(CFM_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin CACHE STRING "Where do you want the output binaries to install to?")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


#Common Library
add_subdirectory(cfm-code)

#Executable for generating fragmentation graphs
include_directories(cfm-code)

#Spectrum Predicter
add_subdirectory(cfm-predict)

add_subdirectory(cfm-boost-test)

add_subdirectory(cfm-train)
