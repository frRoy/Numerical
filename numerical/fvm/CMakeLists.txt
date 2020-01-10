# Numerical fvm
# Library used to solve a pde problem using the finite volume method.
# Copyright 2020, Francois Roy
#=============================================================================
set(HEADERS
  fvm.hpp
  FVMesh.hpp
  FVProblem.hpp
  FVSolver.hpp
  PARENT_SCOPE)

set(SOURCES
  FVMesh.cpp
  FVProblem.cpp
  FVSolver.cpp
  PARENT_SCOPE)

file(GLOB SOURCES *.cpp)
add_library(lfvm SHARED ${SOURCES})

# Setting this to the library only avoids errors in 3rd party dependencies 
# which are compile with 
# e.g. -Werror target_compile_options(example PUBLIC -Wall -pedantic)

# Include code-coverage settings:
target_link_libraries(lfvm PUBLIC coverage_config)