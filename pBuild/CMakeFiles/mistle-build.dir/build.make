# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ynowatzk/repos/mistle

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ynowatzk/repos/mistle/pBuild

# Include any dependencies generated for this target.
include CMakeFiles/mistle-build.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mistle-build.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mistle-build.dir/flags.make

CMakeFiles/mistle-build.dir/src/build_index.cpp.o: CMakeFiles/mistle-build.dir/flags.make
CMakeFiles/mistle-build.dir/src/build_index.cpp.o: ../src/build_index.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ynowatzk/repos/mistle/pBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mistle-build.dir/src/build_index.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mistle-build.dir/src/build_index.cpp.o -c /home/ynowatzk/repos/mistle/src/build_index.cpp

CMakeFiles/mistle-build.dir/src/build_index.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mistle-build.dir/src/build_index.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ynowatzk/repos/mistle/src/build_index.cpp > CMakeFiles/mistle-build.dir/src/build_index.cpp.i

CMakeFiles/mistle-build.dir/src/build_index.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mistle-build.dir/src/build_index.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ynowatzk/repos/mistle/src/build_index.cpp -o CMakeFiles/mistle-build.dir/src/build_index.cpp.s

CMakeFiles/mistle-build.dir/src/indexing_manager.cpp.o: CMakeFiles/mistle-build.dir/flags.make
CMakeFiles/mistle-build.dir/src/indexing_manager.cpp.o: ../src/indexing_manager.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ynowatzk/repos/mistle/pBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/mistle-build.dir/src/indexing_manager.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mistle-build.dir/src/indexing_manager.cpp.o -c /home/ynowatzk/repos/mistle/src/indexing_manager.cpp

CMakeFiles/mistle-build.dir/src/indexing_manager.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mistle-build.dir/src/indexing_manager.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ynowatzk/repos/mistle/src/indexing_manager.cpp > CMakeFiles/mistle-build.dir/src/indexing_manager.cpp.i

CMakeFiles/mistle-build.dir/src/indexing_manager.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mistle-build.dir/src/indexing_manager.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ynowatzk/repos/mistle/src/indexing_manager.cpp -o CMakeFiles/mistle-build.dir/src/indexing_manager.cpp.s

CMakeFiles/mistle-build.dir/src/msp_reader.cpp.o: CMakeFiles/mistle-build.dir/flags.make
CMakeFiles/mistle-build.dir/src/msp_reader.cpp.o: ../src/msp_reader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ynowatzk/repos/mistle/pBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/mistle-build.dir/src/msp_reader.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mistle-build.dir/src/msp_reader.cpp.o -c /home/ynowatzk/repos/mistle/src/msp_reader.cpp

CMakeFiles/mistle-build.dir/src/msp_reader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mistle-build.dir/src/msp_reader.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ynowatzk/repos/mistle/src/msp_reader.cpp > CMakeFiles/mistle-build.dir/src/msp_reader.cpp.i

CMakeFiles/mistle-build.dir/src/msp_reader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mistle-build.dir/src/msp_reader.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ynowatzk/repos/mistle/src/msp_reader.cpp -o CMakeFiles/mistle-build.dir/src/msp_reader.cpp.s

CMakeFiles/mistle-build.dir/src/mgf_reader.cpp.o: CMakeFiles/mistle-build.dir/flags.make
CMakeFiles/mistle-build.dir/src/mgf_reader.cpp.o: ../src/mgf_reader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ynowatzk/repos/mistle/pBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/mistle-build.dir/src/mgf_reader.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mistle-build.dir/src/mgf_reader.cpp.o -c /home/ynowatzk/repos/mistle/src/mgf_reader.cpp

CMakeFiles/mistle-build.dir/src/mgf_reader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mistle-build.dir/src/mgf_reader.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ynowatzk/repos/mistle/src/mgf_reader.cpp > CMakeFiles/mistle-build.dir/src/mgf_reader.cpp.i

CMakeFiles/mistle-build.dir/src/mgf_reader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mistle-build.dir/src/mgf_reader.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ynowatzk/repos/mistle/src/mgf_reader.cpp -o CMakeFiles/mistle-build.dir/src/mgf_reader.cpp.s

CMakeFiles/mistle-build.dir/src/spectrum.cpp.o: CMakeFiles/mistle-build.dir/flags.make
CMakeFiles/mistle-build.dir/src/spectrum.cpp.o: ../src/spectrum.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ynowatzk/repos/mistle/pBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/mistle-build.dir/src/spectrum.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mistle-build.dir/src/spectrum.cpp.o -c /home/ynowatzk/repos/mistle/src/spectrum.cpp

CMakeFiles/mistle-build.dir/src/spectrum.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mistle-build.dir/src/spectrum.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ynowatzk/repos/mistle/src/spectrum.cpp > CMakeFiles/mistle-build.dir/src/spectrum.cpp.i

CMakeFiles/mistle-build.dir/src/spectrum.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mistle-build.dir/src/spectrum.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ynowatzk/repos/mistle/src/spectrum.cpp -o CMakeFiles/mistle-build.dir/src/spectrum.cpp.s

CMakeFiles/mistle-build.dir/src/index_file_writer.cpp.o: CMakeFiles/mistle-build.dir/flags.make
CMakeFiles/mistle-build.dir/src/index_file_writer.cpp.o: ../src/index_file_writer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ynowatzk/repos/mistle/pBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/mistle-build.dir/src/index_file_writer.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mistle-build.dir/src/index_file_writer.cpp.o -c /home/ynowatzk/repos/mistle/src/index_file_writer.cpp

CMakeFiles/mistle-build.dir/src/index_file_writer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mistle-build.dir/src/index_file_writer.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ynowatzk/repos/mistle/src/index_file_writer.cpp > CMakeFiles/mistle-build.dir/src/index_file_writer.cpp.i

CMakeFiles/mistle-build.dir/src/index_file_writer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mistle-build.dir/src/index_file_writer.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ynowatzk/repos/mistle/src/index_file_writer.cpp -o CMakeFiles/mistle-build.dir/src/index_file_writer.cpp.s

CMakeFiles/mistle-build.dir/src/precursor_index.cpp.o: CMakeFiles/mistle-build.dir/flags.make
CMakeFiles/mistle-build.dir/src/precursor_index.cpp.o: ../src/precursor_index.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ynowatzk/repos/mistle/pBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/mistle-build.dir/src/precursor_index.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mistle-build.dir/src/precursor_index.cpp.o -c /home/ynowatzk/repos/mistle/src/precursor_index.cpp

CMakeFiles/mistle-build.dir/src/precursor_index.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mistle-build.dir/src/precursor_index.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ynowatzk/repos/mistle/src/precursor_index.cpp > CMakeFiles/mistle-build.dir/src/precursor_index.cpp.i

CMakeFiles/mistle-build.dir/src/precursor_index.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mistle-build.dir/src/precursor_index.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ynowatzk/repos/mistle/src/precursor_index.cpp -o CMakeFiles/mistle-build.dir/src/precursor_index.cpp.s

CMakeFiles/mistle-build.dir/src/fragment_ion_index.cpp.o: CMakeFiles/mistle-build.dir/flags.make
CMakeFiles/mistle-build.dir/src/fragment_ion_index.cpp.o: ../src/fragment_ion_index.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ynowatzk/repos/mistle/pBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/mistle-build.dir/src/fragment_ion_index.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mistle-build.dir/src/fragment_ion_index.cpp.o -c /home/ynowatzk/repos/mistle/src/fragment_ion_index.cpp

CMakeFiles/mistle-build.dir/src/fragment_ion_index.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mistle-build.dir/src/fragment_ion_index.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ynowatzk/repos/mistle/src/fragment_ion_index.cpp > CMakeFiles/mistle-build.dir/src/fragment_ion_index.cpp.i

CMakeFiles/mistle-build.dir/src/fragment_ion_index.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mistle-build.dir/src/fragment_ion_index.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ynowatzk/repos/mistle/src/fragment_ion_index.cpp -o CMakeFiles/mistle-build.dir/src/fragment_ion_index.cpp.s

CMakeFiles/mistle-build.dir/src/index_file_reader.cpp.o: CMakeFiles/mistle-build.dir/flags.make
CMakeFiles/mistle-build.dir/src/index_file_reader.cpp.o: ../src/index_file_reader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ynowatzk/repos/mistle/pBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/mistle-build.dir/src/index_file_reader.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mistle-build.dir/src/index_file_reader.cpp.o -c /home/ynowatzk/repos/mistle/src/index_file_reader.cpp

CMakeFiles/mistle-build.dir/src/index_file_reader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mistle-build.dir/src/index_file_reader.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ynowatzk/repos/mistle/src/index_file_reader.cpp > CMakeFiles/mistle-build.dir/src/index_file_reader.cpp.i

CMakeFiles/mistle-build.dir/src/index_file_reader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mistle-build.dir/src/index_file_reader.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ynowatzk/repos/mistle/src/index_file_reader.cpp -o CMakeFiles/mistle-build.dir/src/index_file_reader.cpp.s

CMakeFiles/mistle-build.dir/src/configuration.cpp.o: CMakeFiles/mistle-build.dir/flags.make
CMakeFiles/mistle-build.dir/src/configuration.cpp.o: ../src/configuration.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ynowatzk/repos/mistle/pBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/mistle-build.dir/src/configuration.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mistle-build.dir/src/configuration.cpp.o -c /home/ynowatzk/repos/mistle/src/configuration.cpp

CMakeFiles/mistle-build.dir/src/configuration.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mistle-build.dir/src/configuration.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ynowatzk/repos/mistle/src/configuration.cpp > CMakeFiles/mistle-build.dir/src/configuration.cpp.i

CMakeFiles/mistle-build.dir/src/configuration.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mistle-build.dir/src/configuration.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ynowatzk/repos/mistle/src/configuration.cpp -o CMakeFiles/mistle-build.dir/src/configuration.cpp.s

CMakeFiles/mistle-build.dir/src/thread_pool.cpp.o: CMakeFiles/mistle-build.dir/flags.make
CMakeFiles/mistle-build.dir/src/thread_pool.cpp.o: ../src/thread_pool.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ynowatzk/repos/mistle/pBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/mistle-build.dir/src/thread_pool.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mistle-build.dir/src/thread_pool.cpp.o -c /home/ynowatzk/repos/mistle/src/thread_pool.cpp

CMakeFiles/mistle-build.dir/src/thread_pool.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mistle-build.dir/src/thread_pool.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ynowatzk/repos/mistle/src/thread_pool.cpp > CMakeFiles/mistle-build.dir/src/thread_pool.cpp.i

CMakeFiles/mistle-build.dir/src/thread_pool.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mistle-build.dir/src/thread_pool.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ynowatzk/repos/mistle/src/thread_pool.cpp -o CMakeFiles/mistle-build.dir/src/thread_pool.cpp.s

CMakeFiles/mistle-build.dir/src/settings.cpp.o: CMakeFiles/mistle-build.dir/flags.make
CMakeFiles/mistle-build.dir/src/settings.cpp.o: ../src/settings.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ynowatzk/repos/mistle/pBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/mistle-build.dir/src/settings.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mistle-build.dir/src/settings.cpp.o -c /home/ynowatzk/repos/mistle/src/settings.cpp

CMakeFiles/mistle-build.dir/src/settings.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mistle-build.dir/src/settings.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ynowatzk/repos/mistle/src/settings.cpp > CMakeFiles/mistle-build.dir/src/settings.cpp.i

CMakeFiles/mistle-build.dir/src/settings.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mistle-build.dir/src/settings.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ynowatzk/repos/mistle/src/settings.cpp -o CMakeFiles/mistle-build.dir/src/settings.cpp.s

# Object files for target mistle-build
mistle__build_OBJECTS = \
"CMakeFiles/mistle-build.dir/src/build_index.cpp.o" \
"CMakeFiles/mistle-build.dir/src/indexing_manager.cpp.o" \
"CMakeFiles/mistle-build.dir/src/msp_reader.cpp.o" \
"CMakeFiles/mistle-build.dir/src/mgf_reader.cpp.o" \
"CMakeFiles/mistle-build.dir/src/spectrum.cpp.o" \
"CMakeFiles/mistle-build.dir/src/index_file_writer.cpp.o" \
"CMakeFiles/mistle-build.dir/src/precursor_index.cpp.o" \
"CMakeFiles/mistle-build.dir/src/fragment_ion_index.cpp.o" \
"CMakeFiles/mistle-build.dir/src/index_file_reader.cpp.o" \
"CMakeFiles/mistle-build.dir/src/configuration.cpp.o" \
"CMakeFiles/mistle-build.dir/src/thread_pool.cpp.o" \
"CMakeFiles/mistle-build.dir/src/settings.cpp.o"

# External object files for target mistle-build
mistle__build_EXTERNAL_OBJECTS =

mistle-build: CMakeFiles/mistle-build.dir/src/build_index.cpp.o
mistle-build: CMakeFiles/mistle-build.dir/src/indexing_manager.cpp.o
mistle-build: CMakeFiles/mistle-build.dir/src/msp_reader.cpp.o
mistle-build: CMakeFiles/mistle-build.dir/src/mgf_reader.cpp.o
mistle-build: CMakeFiles/mistle-build.dir/src/spectrum.cpp.o
mistle-build: CMakeFiles/mistle-build.dir/src/index_file_writer.cpp.o
mistle-build: CMakeFiles/mistle-build.dir/src/precursor_index.cpp.o
mistle-build: CMakeFiles/mistle-build.dir/src/fragment_ion_index.cpp.o
mistle-build: CMakeFiles/mistle-build.dir/src/index_file_reader.cpp.o
mistle-build: CMakeFiles/mistle-build.dir/src/configuration.cpp.o
mistle-build: CMakeFiles/mistle-build.dir/src/thread_pool.cpp.o
mistle-build: CMakeFiles/mistle-build.dir/src/settings.cpp.o
mistle-build: CMakeFiles/mistle-build.dir/build.make
mistle-build: CMakeFiles/mistle-build.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ynowatzk/repos/mistle/pBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Linking CXX executable mistle-build"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mistle-build.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mistle-build.dir/build: mistle-build

.PHONY : CMakeFiles/mistle-build.dir/build

CMakeFiles/mistle-build.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mistle-build.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mistle-build.dir/clean

CMakeFiles/mistle-build.dir/depend:
	cd /home/ynowatzk/repos/mistle/pBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ynowatzk/repos/mistle /home/ynowatzk/repos/mistle /home/ynowatzk/repos/mistle/pBuild /home/ynowatzk/repos/mistle/pBuild /home/ynowatzk/repos/mistle/pBuild/CMakeFiles/mistle-build.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mistle-build.dir/depend
