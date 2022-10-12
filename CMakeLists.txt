cmake_minimum_required(VERSION 3.23)
project(information_theory)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(information_theory
        simulation_utils.cpp
        simulation_utils.h
        sw_test.cpp encoding_decoding.cpp encoding_decoding.h npy.hpp)
