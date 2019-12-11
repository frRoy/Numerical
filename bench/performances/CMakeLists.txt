# "performances" Target
# run performances
# Copyright 2020, Francois Roy
#=============================================================================
cmake_minimum_required(VERSION 3.10)
project(performances)

file(GLOB SOURCES *.cpp)
add_executable(performances ${SOURCES})
