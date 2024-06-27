cmake_minimum_required(VERSION 3.16)
cmake_policy(VERSION 3.12)
project(KineticGas LANGUAGES C CXX)

SET(CMAKE_CXX_STANDARD 17)

# Find the module development requirements (requires FindPython from 3.17 or scikit-build-core's built-in backport)
set(Python_FIND_VIRTUALENV ONLY)

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -Wno-unused-parameter -Wno-sign-compare -pthread")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG -Wno-unused-parameter -Wno-unused-const-variable")
endif(NOT MSVC)

file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/cpp/*.cpp)

project(pykingas)
set(TARGET libpykingas)

# Add a library using FindPython's tooling (pybind11 also provides a helper like this)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/cpp/external/pybind11")
pybind11_add_module(${TARGET} ${SOURCES} "${CMAKE_CURRENT_SOURCE_DIR}/cpp/Integration/Integration.cpp")

# The install directory is the output (wheel) directory
install(TARGETS ${TARGET} DESTINATION ${CMAKE_SOURCE_DIR}/pykingas)