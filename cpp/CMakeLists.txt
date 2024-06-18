cmake_policy(VERSION 3.12)
project(KineticGas LANGUAGES C CXX)

cmake_minimum_required(VERSION 3.12)
set(PYBIND11_PYTHON_VERSION "3")
string(ASCII 27 Esc)
message("${Esc}[34mPYBIND11_PYTHON_VERSION is : ${PYBIND11_PYTHON_VERSION}")
message("${Esc}[34mPYTHON_EXECUTABLE is : ${PYTHON_EXECUTABLE}")

set(PYBIND11_ROOT "/usr/local/include/pybind11") # This path must contain the pybind11 root dir.
set(JSON_INCLUDE_DIR "/usr/local/include/json/")

if(NOT MSVC)
  if(APPLE)
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare -std=c++17 -pthread -arch arm64")
  else()
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare -std=c++17 -pthread")
  endif(APPLE)
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG -Wno-unused-parameter -Wno-unused-const-variable")
endif(NOT MSVC)

set(SOURCES KineticGas.cpp KineticGas_mthr.cpp Spherical.cpp MieKinGas.cpp
            Integration/Integration.cpp Factorial.cpp)

set(pyTARGET pykingas)
add_subdirectory(${PYBIND11_ROOT} release)
pybind11_add_module(${pyTARGET} ${SOURCES} bindings.cpp)

find_package(Eigen3)
if (Eigen3_FOUND)
  set(libTARGET kineticgas)
  add_library(${libTARGET} ${SOURCES} KineticGas_properties.h)
  target_include_directories(${libTARGET} ${JSON_INCLUDE_DIR})
  target_compile_definitions(${libTARGET} PUBLIC "-DNOPYTHON")
  target_link_libraries(${libTARGET} Eigen3::Eigen)
else()
  message("${Esc}[34mEigen library not found! Not creating target libkineticgas")
endif()
