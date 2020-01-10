# Numerical "unittests" Target
# Tests the numerical package
# Copyright 2020, Francois Roy
#=============================================================================
cmake_minimum_required(VERSION 3.10)
project(numerical-tests)

add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

# Make tests executable
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fdm/FDMeshTest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fdm/FDProblemTest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fdm/ParametersTest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fdm/SparseSolverTest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fem/FEProblemTest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fem/FEMeshTest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fem/FiniteElementTest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fvm/FVMeshTest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fvm/FVProblemTest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fvm/FVSolverTest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/ExecutionTimerTest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/UtilsTest.cpp
  )

add_executable(tests ${TEST_SOURCES})
target_link_libraries(tests PRIVATE Catch lfdm lfem lfvm)

# Test target
add_test(tests tests)
