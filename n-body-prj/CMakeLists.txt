cmake_minimum_required(VERSION 3.21.0)

project(nbody)

add_subdirectory(src)
add_subdirectory(gui)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)