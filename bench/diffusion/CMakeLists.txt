# "diffusion" Target
# run diffusion benchmarks
# Copyright 2020, Francois Roy
#=============================================================================
cmake_minimum_required(VERSION 3.10)
project(diffusion)

file(GLOB SOURCES *.cpp)
add_executable(diffusion ${SOURCES})
target_link_libraries(lfdm)
