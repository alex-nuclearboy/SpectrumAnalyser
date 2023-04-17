cmake_minimum_required(VERSION 3.5)

# Set the project name
project(SpectrumAnalyser VERSION 2.2)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find ROOT package
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Spectrum)

# Add the include directory to the list of include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Include ROOT headers
include(${ROOT_USE_FILE})

# Add executable file to the project
add_executable(run_analysis 
    src/SpectrumAnalyser.cxx 
    src/SpectrumAnalyserCreateHistos.cxx 
    src/SpectrumAnalyserDrawHistos.cxx 
    src/SpectrumAnalyserWriteHistos.cxx
    src/SpectrumAnalyserPeakFinder.cxx
    src/SpectrumAnalyserPrecalib.cxx
    run_analysis.cpp)

# Link ROOT libraries
target_link_libraries(run_analysis PRIVATE ${ROOT_LIBRARIES})

# Set output directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

# Additional compile options
add_compile_options(run_analysis PRIVATE -Wall -Wextra -pedantic)

