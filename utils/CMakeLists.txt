# Numerical "utils" Target
# Library utils fuctions
# Copyright 2020, Francois Roy
#=============================================================================
set(HEADERS
  ExecutionTimer.hpp
  Utils.hpp
  )

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${HEADERS}
     DESTINATION ${NUMERICAL_INCLUDE_DIR}/utils)

# file(GLOB SOURCES *.cpp)
# add_library(lutils SHARED ${SOURCES})
