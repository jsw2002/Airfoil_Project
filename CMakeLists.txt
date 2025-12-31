# Minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# The name of your project
project(AirfoilCFD CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directories
include_directories(src)
include_directories(tests)

# Main Solver Executable
add_executable(flow_solver 
    src/main.cpp
    src/simulation.cpp
    src/simulation.h
)

# Test Executable
add_executable(run_tests
    tests/main_test.cpp
    tests/test_simulation.cpp
    src/simulation.cpp
    src/simulation.h
)