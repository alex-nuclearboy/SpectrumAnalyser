cmake_minimum_required(VERSION 3.5)

# Set the project name
project(SpectrumAnalyser)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find ROOT package
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Spectrum)

# Set include directories
include_directories(${ROOT_INCLUDE_DIRS})

# Include ROOT headers
include(${ROOT_USE_FILE})

# Add class definition file to the project
add_library(SpectrumAnalyser OBJECT SpectrumAnalyser.cxx SpectrumAnalyser.h)

# Add executable file to the project
add_executable(run_analysis run_analysis.cpp $<TARGET_OBJECTS:SpectrumAnalyser>)

# Link ROOT libraries
target_link_libraries(run_analysis ${ROOT_LIBRARIES})

# Set output directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

# Additional compile options
add_compile_options(-Wall -Wextra -pedantic)

message("ROOT Include Directories: " ${ROOT_INCLUDE_DIRS})
message("ROOT Include Files: " ${ROOT_USE_FILE})
message("ROOT Libraries: " ${ROOT_LIBRARIES})
