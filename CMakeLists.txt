# Numerical
# A code used to showcase numerical methods.
# Copyright 2020, Francois Roy
#=============================================================================
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(Numerical CXX)
set(PROJECT_VERSION_RELEASE 0)
set(PROJECT_VERSION_MAJOR "2019")
set(PROJECT_VERSION_MINOR "0")
set(PROJECT_VERSION_PATCH "0")
set(PROJECT_DESCRIPTION  "A code used to showcase numerical methods.")
set(PROJECT_AUTHOR  "Francois Roy")
set(PROJECT_URL  "https://bitbucket.org/frnsroy/numerical")
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check compiler version

# Check for GCC version - earlier versions have insuffcient C++11
# support, or bugs.
if (CMAKE_COMPILER_IS_GNUCXX)
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8.3)
    message(FATAL_ERROR "GCC version must be at least 4.8.3 
            (for sufficient C++11  support and to avoid some compiler bugs). 
            You have version ${CMAKE_CXX_COMPILER_VERSION}")
  endif()
endif()

# set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)  # for valgrind target
set(PROJECT_ENABLE_CODE_COVERAGE OFF)  # for code coverage
# Debugging Options
set (CMAKE_VERBOSE_MAKEFILE 0) # Use 1 for debugging, 0 otherwise.

set(NUMERICAL_INCLUDE_DIR ${PROJECT_BINARY_DIR}/include)
include_directories(${NUMERICAL_INCLUDE_DIR})

# Add flags for generating code coverage reports
if (PROJECT_ENABLE_CODE_COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage 
      -fPIC -O0")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs 
      -ftest-coverage -fPIC -O0")
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_CONFIG_PATH ${PROJECT_SOURCE_DIR}/cmake.config)

set (platform)
if (APPLE)
  set (platform "apple")
  if(NOT CMAKE_OSX_ARCHITECTURES)
    set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "build for 64-bit")
  endif()
elseif (UNIX)
  set (platform "unix")
elseif (WIN32)
  set (platform "win32")
  message(FATAL_ERROR 
      "This project is not supported on ${platform} platforms."
      )
endif()

find_package(Git REQUIRED)
set(EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/third_party)
# find_package(Eigen3 REQUIRED)
find_package(spdlog REQUIRED)
find_package(Catch2 REQUIRED)

message("-- The Eigen3 include dir is: ${EIGEN3_INCLUDE_DIR}")
message("-- The Cathc2 include dir is: ${Catch2_INCLUDE_DIRS}")
message("-- The spdlog include dir is: ${spdlog_INCLUDE_DIR}")

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Catch2_INCLUDE_DIRS})
include_directories(${spdlog_INCLUDE_DIR})

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${CMAKE_CONFIG_PATH}/NumericalConfig.h.in"
  "${NUMERICAL_INCLUDE_DIR}/NumericalConfig.h"
  )

add_subdirectory(numerical)
add_subdirectory(utils)

# Add documentation
# target doc
# set files manually (for now)
set(INPUT_DOC  
  ${PROJECT_SOURCE_DIR}/doc
  ${PROJECT_SOURCE_DIR}/numerical
  ${PROJECT_SOURCE_DIR}/numerical/fdm
  ${PROJECT_SOURCE_DIR}/utils
   ${PROJECT_SOURCE_DIR}/bench/diffusion
  )
string (REPLACE ";" " " INPUT_DOC "${INPUT_DOC}")
add_subdirectory(doc)

# Add tests and benchmarks

# target diffusion
add_subdirectory(${PROJECT_SOURCE_DIR}/bench/diffusion EXCLUDE_FROM_ALL)

# tests coverage?
enable_testing()
# Add target "tests", but do not add to default target
add_subdirectory(${PROJECT_SOURCE_DIR}/tests EXCLUDE_FROM_ALL)
include(CTest)

#TODO generate properflags for gcovr --> generate .gcda file #################
# target coverage
if (PROJECT_ENABLE_CODE_COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)
  include(CodeCoverage)  # doesn't generate the .gcda file
  APPEND_COVERAGE_COMPILER_FLAGS()
  SETUP_TARGET_FOR_COVERAGE_GCOVR_HTML(
    NAME coverage
    EXECUTABLE ctest
    DEPENDENCIES tests
    )
endif()
##############################################################################

# target valgrind
find_program(VALGRIND "valgrind")
if(VALGRIND)
    add_custom_target(valgrind
        COMMAND "${VALGRIND}" --tool=memcheck --leak-check=yes 
            --show-reachable=yes --num-callers=20 
            --track-fds=yes $<TARGET_FILE:tests>)
endif()

# target performances
add_subdirectory(${PROJECT_SOURCE_DIR}/bench/performances EXCLUDE_FROM_ALL)
