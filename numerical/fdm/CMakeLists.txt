# Numerical "fdm" Target
# Executable used to solve a pde problem using finite difference.
# Copyright 2020, Francois Roy
#=============================================================================
set(HEADERS
  fdm.hpp
  Mesh.hpp
  Parameters.hpp
  Problem.hpp
  SparseSolver.hpp
  PARENT_SCOPE)

set(SOURCES
  Mesh.cpp
  Parameters.cpp
  Problem.cpp
  SparseSolver.cpp
  PARENT_SCOPE)

file(GLOB SOURCES *.cpp)
add_library(lfdm SHARED ${SOURCES})
