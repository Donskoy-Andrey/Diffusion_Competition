cmake_minimum_required(VERSION 2.14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2")

project(Diffusion DESCRIPTION "Sirius competition" LANGUAGES CXX)

find_package(MPI REQUIRED)

include_directories(
source
source/mesh
)

add_executable(
    main source/main.cpp
    source/mesh/creating.h
)

target_link_libraries(main MPI::MPI_CXX)