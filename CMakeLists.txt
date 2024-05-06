cmake_minimum_required(VERSION 3.16)
cmake_policy(VERSION 3.12)
project(KineticGas LANGUAGES C CXX)

SET(CMAKE_CXX_STANDARD 17)

if(NOT MSVC)
  # set(CMAKE_CXX_COMPILER /usr/bin/g++) # Setting these can, on some systems, cause cmake to go into an infinite loop (see: build.sh to set the compiler)
  # set(CMAKE_C_COMPILER /usr/bin/gcc) # Setting these can, on some systems, cause cmake to go into an infinite loop (see: build.sh to set the compiler)
  if(APPLE)
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare -pthread -arch arm64")
  else()
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare -pthread")
  endif(APPLE)
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG -Wno-unused-parameter -Wno-unused-const-variable")
endif(NOT MSVC)

set(SOURCES cpp/KineticGas.cpp cpp/KineticGas_mthr.cpp cpp/Spherical.cpp cpp/MieKinGas.cpp
            cpp/Integration/Integration.cpp cpp/Factorial.cpp cpp/bindings.cpp)

if (CMAKE_BUILD_TYPE MATCHES Release)
    project(pykingas)
    set(TARGET pykingas)
else()
    project(pykingas_d)
    set(TARGET pykingas_d)
endif()

# Find the module development requirements (requires FindPython from 3.17 or
# scikit-build-core's built-in backport)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Add a library using FindPython's tooling (pybind11 also provides a helper like
# this)
python_add_library(${TARGET} MODULE ${SOURCES} WITH_SOABI)
target_link_libraries(${TARGET} PRIVATE pybind11::headers)

# The install directory is the output (wheel) directory
install(TARGETS ${TARGET} DESTINATION pykingas)