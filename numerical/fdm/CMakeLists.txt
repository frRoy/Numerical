# Numerical fdm
# Library used to solve a pde problem using the finite difference method.
# Copyright 2020, Francois Roy
#=============================================================================
set(HEADERS
  fdm.hpp
  FDMesh.hpp
  Parameters.hpp
  FDProblem.hpp
  SparseSolver.hpp
  PARENT_SCOPE)

set(SOURCES
  FDMesh.cpp
  Parameters.cpp
  FDProblem.cpp
  SparseSolver.cpp
  PARENT_SCOPE)

file(GLOB SOURCES *.cpp)
add_library(lfdm SHARED ${SOURCES})

# Setting this to the library only avoids errors in 3rd party dependencies 
# which are compile with 
# e.g. -Werror target_compile_options(example PUBLIC -Wall -pedantic)

# Include code-coverage settings:
target_link_libraries(lfdm PUBLIC coverage_config)
