# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.5)

# Define the project and its version
project(SpectrumAnalyser)

# Add the include directory to the list of include directories
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)
include_directories(${INCLUDE_DIR})

# Add the source files to the project
add_library(SpectrumAnalyser 
    SpectrumAnalyser.cxx
    SpectrumAnalyserCreateHistos.cxx
    SpectrumAnalyserDrawHistos.cxx
    SpectrumAnalyserWriteHistos.cxx)

include_directories(${INCLUDE_DIR})

# Link the ROOT libraries to the executable
target_link_libraries(SpectrumAnalyser PUBLIC ${ROOT_LIBRARIES})
