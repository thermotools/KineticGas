cmake_minimum_required(VERSION 3.16)
cmake_policy(VERSION 3.12)
project(KineticGas LANGUAGES C CXX)

string(ASCII 27 Esc)
set(ColorDefault "${Esc}[0m")
set(ColorBlue "${Esc}[34m")
set(ColorYellow "${Esc}[93m")
set(ColorRed "${Esc}[31m")
set(ColorDefault "${Esc}[0m")

option(build "Recompile KineticGas?" ON) # If set to off, assumes that KineticGas is installed and does nothing. Useful in combination with pip
if(NOT build)
  message(STATUS "Assuming KineticGas is already built and installed. Exiting ...")
  message(STATUS "${ColorBlue}To disable this behaviour, clear the CMakeCache and re-run cmake.${ColorDefault}")
  return()
endif()

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
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG -Wno-unused-parameter -Wno-unused-const-variable")
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

if (DEFINED ENV{THERMOPACK_DIR})
  set(THERMOPACK_DIR $ENV{THERMOPACK_DIR})
  message(STATUS "Searching for thermopack in THERMOPACK_DIR: $ENV{THERMOPACK_DIR}")
else()
  message(STATUS "THERMOPACK_DIR is not defined!")
  set(THERMOPACK_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp/external/thermopack)
  message(STATUS "Setting THERMOPACK_DIR: ${THERMOPACK_DIR}")
endif()

find_package(THERMOPACK REQUIRED)
if(NOT THERMOPACK_INSTALLED)
  message(STATUS "${ColorRed}ThermoPack was found, but is not installed. If you have installed ThermoPack, 
set the environment variable 'THERMOPACK_DIR=</path/to/thermopack>' to help cmake find your installation.${ColorDefault}")
  message(STATUS "${ColorBlue}Target 'thermopack' created from the files found at ${THERMOPACK_ROOT}.${ColorDefault}")
  message(STATUS "${ColorBlue}Try Running `make install` to build and install thermopack, before re-running `cmake .. && make install` to build KineticGas.${ColorDefault}")
  add_subdirectory(${THERMOPACK_ROOT}) # We need to build thermopack, because it has not been installed.
  return() # Hoping that thermopack target has been created by thermopack, exiting ...
endif()

option(purecpp "Build pure cpp lib" OFF)
if (purecpp)
  message(STATUS "${ColorBlue}Generating C++ library targets. Run `cmake [-Dpurecpp=OFF] ..` to skip this target (you will need to clear the CMakeCache first).${ColorDefault}")
  
  if(NOT DEFINED(FLUID_DIR))
    set(FLUID_DIR "../fluids" CACHE STRING "Default search path for fluid files")
  endif()
  message(STATUS "Default search path for fluid files (relative to dynamic library): ${FLUID_DIR}")

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
  message(STATUS "${ColorYellow}Not generating C++ library targets. Run `cmake -Dpurecpp=ON ..` to generate these targets (you will need to clear the CMakeCache first).${ColorDefault}")
endif()

option(pylib "Build library for pykingas" ON)
if (pylib)
  message(STATUS "${ColorBlue}Generating libpykingas target. Run `cmake -Dpylib=OFF ..` to skip this target (you will need to clear the CMakeCache first).${ColorDefault}")
  set(pyTARGET libpykingas)
  set(pykingas_path ${CMAKE_CURRENT_SOURCE_DIR}/pykingas)
  add_subdirectory(${PYBIND11_ROOT})
  pybind11_add_module(${pyTARGET} ${COMMON_SRC} ${SRC_DIR}/bindings.cpp)

  if(NOT MSVC)
    set_target_properties(${pyTARGET} PROPERTIES 
                          INSTALL_RPATH ${THERMOPACK_ROOT}/installed
                          INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif()

  set(pyFLUID_DIR "./fluids" CACHE STRING "Default search path for fluid files (pykingas)")
  message(STATUS "pykingas search path for fluids (relative to pykingas directory): ${pyFLUID_DIR}")

  target_include_directories(${pyTARGET} PUBLIC ${JSON_INCLUDE} ${EIGEN_INCLUDE} ${AUTODIFF_INCLUDE})
  target_link_libraries(${pyTARGET} PRIVATE thermopack)
  target_compile_definitions(${pyTARGET} PUBLIC "-DPYLIB -DFLUID_DIR=${pyFLUID_DIR}")

  install(TARGETS ${pyTARGET} DESTINATION ${pykingas_path})
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/fluids DESTINATION ${pykingas_path})
  # file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/fluids DESTINATION ${pykingas_path})
endif()