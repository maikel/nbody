cmake_minimum_required(VERSION 3.13)
project(nbody C CXX)

# Set a default build type if none was specified
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(default_build_type "Debug")
endif()
 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
else()
  message(STATUS "Setting build type to '${default_build_type}'.")
endif()

# Detect if we use conan to build this project and its dependencies.
# If yes, we include the neccessary build information, i.e. include paths and
# linker options.
if (EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  conan_basic_setup()
endif()

find_package(Boost REQUIRED COMPONENTS program_options container)
find_package(Vc 1.4 REQUIRED)
find_package(fmt REQUIRED)

add_executable(n-body main.cpp run_simulation.cpp)
target_link_libraries(n-body Boost::program_options Boost::container Vc::Vc fmt::fmt)
target_compile_features(n-body PRIVATE cxx_std_14)
