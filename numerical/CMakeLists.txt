# Numerical
# Set sources and headers
# Copyright 2020, Francois Roy
#=============================================================================
set(NUMERICAL_H numerical.hpp)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${NUMERICAL_H}
     DESTINATION ${NUMERICAL_INCLUDE_DIR}/numerical)

# All files and directories in this directory
file(GLOB NUMERICAL_CONTENT RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *)

# Initialize list
set(NUMERICAL_DIRS)

# List of excluded files and directories
set(EXCLUDES ${NUMERICAL_H} CMakeLists.txt)

# Iterate over all directories and files and append to NUMERICAL_DIRS
foreach(_DIR_NAME ${NUMERICAL_CONTENT})
  list(FIND EXCLUDES ${_DIR_NAME} INDEX_EXCLUDES)
  if (${INDEX_EXCLUDES} LESS 0)
    list(APPEND NUMERICAL_DIRS ${_DIR_NAME})
  endif()
endforeach()

# Initialize lists
set(NUMERICAL_HEADERS)
set(NUMERICAL_SOURCES)

foreach(DIR ${NUMERICAL_DIRS})
  # Each subdirectory defines HEADERS and SOURCES
  add_subdirectory(${DIR})
  set(HEADERS_FULL_PATH)
  foreach(HEADER_FILE ${HEADERS})
    list(APPEND HEADERS_FULL_PATH 
         ${CMAKE_CURRENT_SOURCE_DIR}/${DIR}/${HEADER_FILE})
  endforeach()
  file(COPY ${HEADERS_FULL_PATH} DESTINATION 
       ${NUMERICAL_INCLUDE_DIR}/numerical/${DIR})
  list(APPEND NUMERICAL_HEADERS ${HEADERS_FULL_PATH})

  foreach(SOURCE_FILE ${SOURCES})
    list(APPEND NUMERICAL_SOURCES 
         ${CMAKE_CURRENT_SOURCE_DIR}/${DIR}/${SOURCE_FILE})
  endforeach()
endforeach()


