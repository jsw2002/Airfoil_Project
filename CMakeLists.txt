# Minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# The name of your project (e.g., AirfoilCFD)
project(AirfoilCFD CXX)

# Set the C++ standard (modern C++ is recommended)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Define your main executable:
# "airfoil_solver" is the name of the program it will build
# "main.cpp" is the source file(s) needed to build it
add_executable(airfoil_solver main.cpp)