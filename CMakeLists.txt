cmake_minimum_required(VERSION 3.16)
cmake_policy(VERSION 3.12)
project(KineticGas LANGUAGES C CXX)

string(ASCII 27 Esc)
set(ColorBlue "${Esc}[34m")
set(ColorRed "${Esc}[31m")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "No build type specified, defaulting to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

SET(CMAKE_CXX_STANDARD 17)

# Find the module development requirements (requires FindPython from 3.17 or scikit-build-core's built-in backport)
set(Python_FIND_VIRTUALENV ONLY)

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "-Wfatal-errors -Wall -Wextra -Wno-unknown-pragmas -Wno-unused-parameter -Wno-sign-compare -pthread")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -D_GLIBCXX_DEBUG -Wno-unused-parameter -Wno-unused-const-variable")
endif(NOT MSVC)

list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cpp/external)

set(JSON_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/external/json/include/)
set(EIGEN_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/external/eigen/)
set(AUTODIFF_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/external/autodiff)
set(PYBIND11_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/cpp/external/pybind11)

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp)
file(GLOB COMMON_SRC ${SRC_DIR}/*.cpp)
list(REMOVE_ITEM COMMON_SRC ${SRC_DIR}/bindings.cpp ${SRC_DIR}/run_kineticgas.cpp)
list(APPEND COMMON_SRC ${SRC_DIR}/Integration/Integration.cpp)

find_package(thermopack REQUIRED)
if(NOT THERMOPACK_INSTALLED)
  message("${ColorRed}ThermoPack was found, but is not installed. Resorting to building ThermoPack from source. If you have installed ThermoPack, 
          set the environment variable 'THERMOPACK_DIR=</path/to/thermopack>' to help cmake find your installation.")
  add_subdirectory(${THERMOPACK_ROOT}) # We need to build thermopack, because it has not been installed.
endif()

set(FLUID_DIR ${CMAKE_CURRENT_SOURCE_DIR}/fluids CACHE STRING "Default search path for fluid files")
message(STATUS "Default search path for fluid files : ${FLUID_DIR}")

option(purecpp "Build pure cpp lib" OFF)
if (purecpp)
  set(libTARGET kineticgas)
  add_library(${libTARGET} SHARED ${COMMON_SRC})
  target_include_directories(${libTARGET} PUBLIC ${JSON_INCLUDE} ${EIGEN_INCLUDE} ${AUTODIFF_INCLUDE})
  target_compile_definitions(${libTARGET} PUBLIC "-DFLUID_DIR=${FLUID_DIR}")
  target_link_libraries(${libTARGET} thermopack)
  install(TARGETS ${libTARGET} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)

  set(runTARGET run_kineticgas)
  add_executable(${runTARGET} ${SRC_DIR}/run_kineticgas.cpp)
  target_link_libraries(${runTARGET} ${libTARGET})
else()
  message("${ColorBlue}Not generating C++ library targets. Run `cmake -Dpurecpp=ON ..` to generate these targets (you will need to clear the CMakeCache first).")
endif()

option(pylib "Build library for pykingas" ON)
if (pylib)
  message("${ColorBlue}Generating libpykingas target. Run `cmake -Dpylib=OFF ..` to skip this target (you will need to clear the CMakeCache first).")
  set(pyTARGET libpykingas)
  add_subdirectory(${PYBIND11_ROOT})
  pybind11_add_module(${pyTARGET} ${COMMON_SRC} ${SRC_DIR}/bindings.cpp)
  target_include_directories(${pyTARGET} PUBLIC ${JSON_INCLUDE} ${EIGEN_INCLUDE} ${AUTODIFF_INCLUDE})
  target_link_libraries(${pyTARGET} PRIVATE thermopack)
  target_compile_definitions(${pyTARGET} PUBLIC "-DPYLIB -DFLUID_DIR=${FLUID_DIR}")

  install(TARGETS ${pyTARGET} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/pykingas)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/fluids DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/pykingas/)
endif()