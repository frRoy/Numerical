# Numerical fem
# Library used to solve a pde problem using finite element.
# Copyright 2020, Francois Roy
#=============================================================================
set(HEADERS
  fem.hpp
  FEProblem.hpp
  FiniteElement.hpp
  FEMesh.hpp
  PARENT_SCOPE)

set(SOURCES
  FEProblem.cpp
  FiniteElement.cpp
  FEMesh.cpp
  PARENT_SCOPE)

file(GLOB SOURCES *.cpp)
add_library(lfem SHARED ${SOURCES})

# Setting this to the library only avoids errors in 3rd party dependencies 
# which are compile with 
# e.g. -Werror target_compile_options(example PUBLIC -Wall -pedantic)

# Include code-coverage settings:
target_link_libraries(lfem PUBLIC coverage_config)