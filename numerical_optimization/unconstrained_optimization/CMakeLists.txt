cmake_minimum_required(VERSION 3.10)

project(cubic_spline_smooth)

set(CMAKE_BUILD_TYPE "Debug")

include_directories(/mnt/d/Github/AutoDrive/eigen)

add_executable(cubic_smooth main.cpp)