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

find_package(Boost 1.86 REQUIRED COMPONENTS ${Boost_components})
include_directories(${Boost_INCLUDE_DIRS})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


include(ExternalProject)
set(EXTERNAL_INSTALL_LOCATION ${CMAKE_SOURCE_DIR}/rdkit)
ExternalProject_Add(RDKit
		PREFIX ${EXTERNAL_INSTALL_LOCATION}
		GIT_REPOSITORY https://github.com/rdkit/rdkit
		CMAKE_ARGS -DRDK_BUILD_PYTHON_WRAPPERS=OFF -DRDK_BUILD_INCHI_SUPPORT=ON
)
include_directories(${EXTERNAL_INSTALL_LOCATION}/src/RDKit/Code)
include_directories(${EXTERNAL_INSTALL_LOCATION}/src/RDKit/External)
link_directories(${EXTERNAL_INSTALL_LOCATION}/src/RDKit-build/lib)


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
link_directories(LPSOLVE_LIBRARY)
mark_as_advanced(LPSOLVE_LIBRARY)

find_package(Threads)
find_package(OpenMP)

#Common Library
add_subdirectory(cfm-code)
