# Numerical "pynumerical" Target
# Python binding for the numerical package
# Copyright 2020, Francois Roy
#=============================================================================
cmake_minimum_required(VERSION 3.10)
project(pynumerical)

#TODO move to setup.py
# find_package(pybind11 REQUIRED)
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/pybind11)
set(PYBIND11_CPP_STANDARD -std=c++14)
# message("-- The pybin11 directory is: ${pybind11_DIR}")

set(PYTHON_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
  )

pybind11_add_module(pynumerical ${PYTHON_SOURCES})
