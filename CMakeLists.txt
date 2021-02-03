### example CMakeLists.txt to develop C++ programs using OpenMS
project(simple_example)
cmake_minimum_required(VERSION 3.13)  # CMake version check
set(CMAKE_CXX_STANDARD 20)            # Enable c++20 standard


## list all your executables here (a corresponding .cpp file should exist, e.g. Main.cpp)
set(my_executables
        main
        )

## list all classes here, which are required by your executables
## (all these classes will be linked into a library)
set(my_sources
        src/main.cpp src/spectrum.cpp src/msp_reader.cpp src/scores.cpp src/library.cpp src/library.h src/mgf_reader.cpp src/mgf_reader.h src/spectral_search.cpp src/spectral_search.h src/match.cpp src/match.h src/fragment_ion_index.cpp src/fragment_ion_index.h src/precursor_index.cpp src/precursor_index.h
        )

## find OpenMS configuration and register target "OpenMS" (our library)
#find_package(OpenMS)
## if the above fails you can try calling cmake with -D OpenMS_DIR=/path/to/OpenMS/
## or modify the find_package() call accordingly
find_package(OpenMS PATHS "/home/ynowatzk/repos/OpenMS-build")
#find_package(qt5 PATHS "/usr/lib/x86_64-linux-gnu/qt5")
# TODO these had to be added, for Qt5 and OpenMP to be found
find_package(Qt5 COMPONENTS Core Network REQUIRED)
#find_package(Qt5 REQUIRED)
find_package(OpenMP)

# check whether the OpenMS package was found
if (OpenMS_FOUND)

    # check if the variable containing the include directories is defined
    if(NOT OpenMS_INCLUDE_DIRECTORIES)
        set(_message "The variable \${OpenMS_INCLUDE_DIRECTORIES} is not defined.")
        set(_message "${_message}This CMakeLists.txt and your build of OpenMS seem incompatible.")
        set(_message "${_message}Please use the latest version from the OpenMS release!")
        message(FATAL_ERROR ${_message})
    endif()

    ## include directories for OpenMS headers (and contrib)
    # Note: If you want to link against a specific library contained
    #       in the OpenMS package you should also list the
    #       corresponding include variable here, e.g.,
    #       OpenMS_GUI -> ${OpenMS_GUI_INCLUDE_DIRECTORIES}
    include_directories(${OpenMS_INCLUDE_DIRECTORIES})

    ## append precompiler macros and compiler flags specific to OpenMS
    ## Warning: this could be harmful to your project. Check this if problems occur.
    ## Also, use this to add your own compiler flags, e.g. for OpenMP support.
    ## e.g. for Visual Studio use /openmp
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMS_ADDCXX_FLAGS} -fPIC -O3") #todo added -fPIC and -O3
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMS_ADDCXX_FLAGS} /openmp")

    add_definitions(/DBOOST_ALL_NO_LIB) ## disable auto-linking of boost libs (Boost tends to guess wrong lib names)

    ## library with additional classes from above
    add_library(my_custom_lib STATIC ${my_sources})

    ## add targets for the executables
    foreach(i ${my_executables})
        add_executable(${i} src/${i}.cpp)
        ## link executables against OpenMS
        target_link_libraries(${i} OpenMS my_custom_lib)
    endforeach(i)


else(OpenMS_FOUND)
    message(FATAL_ERROR "OpenMSConfig.cmake file not found!")
endif(OpenMS_FOUND)