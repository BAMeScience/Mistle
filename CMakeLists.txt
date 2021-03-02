cmake_minimum_required(VERSION 3.13)  # CMake version check
project(simple_example)               # Create project "simple_example"
set(CMAKE_CXX_STANDARD 20)            # Enable c++20 standard
set(CMAKE_CXX_FLAGS "-O3")

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES src/spectrum.cpp src/msp_reader.cpp src/scores.cpp src/library.cpp src/library.h src/mgf_reader.cpp src/mgf_reader.h src/spectral_search.h src/match.cpp src/match.h src/fragment_ion_index.cpp src/fragment_ion_index.h src/precursor_index.cpp src/precursor_index.h src/index_file_reader.h src/index_file_reader.cpp src/index_file_writer.cpp src/index_file_writer.h) #deleted src/spectral_search.cpp

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(test_executable src/main.cpp ${SOURCE_FILES})
add_executable(quick_scan src/quick_scan.cpp ${SOURCE_FILES} src/scanner.cpp src/scanner.h)
add_executable(mistle-build src/build_index.cpp src/indexing_manager.cpp src/indexing_manager.h src/msp_reader.cpp src/msp_reader.h src/spectrum.h src/spectrum.cpp src/index_file_writer.cpp src/index_file_writer.h src/precursor_index.cpp src/precursor_index.h src/fragment_ion_index.cpp src/fragment_ion_index.h src/index_file_reader.cpp src/index_file_reader.h)