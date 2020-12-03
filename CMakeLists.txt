cmake_minimum_required(VERSION 3.13)  # CMake version check
project(simple_example)               # Create project "simple_example"
set(CMAKE_CXX_STANDARD 20)            # Enable c++20 standard

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES src/main.cpp src/spectrum.cpp src/msp_reader.cpp src/scores.cpp src/library.cpp src/library.h src/mgf_reader.cpp src/mgf_reader.h)

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(test_executable ${SOURCE_FILES})