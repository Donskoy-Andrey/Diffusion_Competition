cmake_minimum_required(VERSION 2.14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

project(Diffusion DESCRIPTION "Sirius competition" LANGUAGES CXX)

include_directories(
source
source/mesh
)

add_executable(
    main source/main.cpp
    source/mesh/creating.h
)