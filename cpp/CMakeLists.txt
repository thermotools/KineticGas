project(KineticGas LANGUAGES C CXX)

cmake_minimum_required(VERSION 3.12)
set(PYBIND11_PYTHON_VERSION "3.11")
string(ASCII 27 Esc)
message("${Esc}[34mPYBIND11_PYTHON_VERSION is : ${PYBIND11_PYTHON_VERSION}")
message("${Esc}[34mPYTHON_EXECUTABLE is : ${PYTHON_EXECUTABLE}")

set(PYBIND11_ROOT "/usr/local/include/pybind11") # This path must contain the pybind11 root dir.

if(NOT MSVC)
  # set(CMAKE_CXX_COMPILER /usr/bin/g++) # Setting these can, on some systems, cause cmake to go into an infinite loop (see: build.sh to set the compiler)
  # set(CMAKE_C_COMPILER /usr/bin/gcc) # Setting these can, on some systems, cause cmake to go into an infinite loop (see: build.sh to set the compiler)
  if(APPLE)
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare -std=c++17 -pthread -arch arm64")
  else()
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare -std=c++17 -pthread")
  endif(APPLE)
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG -Wno-unused-parameter -Wno-unused-const-variable")
endif(NOT MSVC)

set(SOURCES KineticGas.cpp KineticGas_mthr.cpp Spherical.cpp MieKinGas.cpp Sutherland.cpp QuantumMie.cpp
            Integration/Integration.cpp Factorial.cpp bindings.cpp)

if (CMAKE_BUILD_TYPE MATCHES Release)
    project(KineticGas_r)
    set(TARGET KineticGas_r)
else()
    project(KineticGas_d)
    set(TARGET KineticGas_d)
endif()

add_subdirectory(${PYBIND11_ROOT} release)
pybind11_add_module(${TARGET} ${SOURCES})