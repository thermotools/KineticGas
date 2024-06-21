cmake_policy(VERSION 3.12)
project(KineticGas LANGUAGES C CXX)

cmake_minimum_required(VERSION 3.12)
set(PYBIND11_PYTHON_VERSION "3")
string(ASCII 27 Esc)
message("${Esc}[34mPYBIND11_PYTHON_VERSION is : ${PYBIND11_PYTHON_VERSION}")
message("${Esc}[34mPYTHON_EXECUTABLE is : ${PYTHON_EXECUTABLE}")

set(PYBIND11_ROOT "/usr/local/include/pybind11") # This path must contain the pybind11 root dir (https://github.com/pybind/pybind11)
set(JSON_INCLUDE "/usr/local/include/") # The directory json/include/nlohmann/ is here (https://github.com/nlohmann/json)
set(THERMOPACK_LIB "/usr/local/lib/libthermopack_optim_gfortran_Darwin.dylib") # This is where thermopack is
set(THERMOPACK_INCLUDE "/usr/local/include/") # The cppThermopack directory is here
set(FLUID_DIR "${CMAKE_SOURCE_DIR}/../fluids") # Path to fluid file directory

if(NOT MSVC)
  if(APPLE)
    set(CMAKE_CXX_FLAGS "-Wfatal-errors -Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare -std=c++17 -pthread")
  else()
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare -std=c++17 -pthread")
  endif(APPLE)
  set(CMAKE_CXX_FLAGS_RELEASE "-g -O0 -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG -Wno-unused-parameter -Wno-unused-const-variable")
endif(NOT MSVC)

set(SOURCES KineticGas.cpp KineticGas_mthr.cpp Spherical.cpp MieKinGas.cpp
            Integration/Integration.cpp Factorial.cpp)

set(pyTARGET libpykingas)
add_subdirectory(${PYBIND11_ROOT} release)
pybind11_add_module(${pyTARGET} ${SOURCES} bindings.cpp)

find_package(Eigen3)
if (Eigen3_FOUND)
  set(libTARGET kineticgas)
  add_library(${libTARGET} STATIC ${SOURCES} KineticGas_properties.cpp global_params.cpp)
  target_include_directories(${libTARGET} PUBLIC ${JSON_INCLUDE_DIR} ${THERMOPACK_INCLUDE})
  target_compile_definitions(${libTARGET} PUBLIC "-DNOPYTHON -DFLUID_DIR=${FLUID_DIR}")
  target_link_libraries(${libTARGET} Eigen3::Eigen ${THERMOPACK_LIB})

  set(runTARGET run_kineticgas)
  add_executable(${runTARGET} run_kineticgas.cpp)
  target_link_libraries(${runTARGET} ${libTARGET})
else()
  message("${Esc}[34mEigen library not found! Not creating targets libkineticgas and run_kineticgas")
endif()
