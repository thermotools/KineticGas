cmake_minimum_required(VERSION 3.16)
cmake_policy(VERSION 3.12)
project(KineticGas LANGUAGES C CXX)

SET(CMAKE_CXX_STANDARD 17)

# Find the module development requirements (requires FindPython from 3.17 or scikit-build-core's built-in backport)
set(Python_FIND_VIRTUALENV ONLY) 
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
message("Python is : ${Python_EXECUTABLE}")
find_package(pybind11 CONFIG)
if (pybind11::FOUND)
  set(pybind11_HEADERS ${pybind11::headers})
else()
  set(pybind11_HEADERS "/usr/local/include/pybind11/include")
endif()

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -Wno-unused-parameter -Wno-sign-compare -pthread")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG -Wno-unused-parameter -Wno-unused-const-variable")
endif(NOT MSVC)

set(SOURCES cpp/KineticGas.cpp cpp/KineticGas_mthr.cpp cpp/Spherical.cpp cpp/MieKinGas.cpp
            cpp/Integration/Integration.cpp cpp/Factorial.cpp cpp/bindings.cpp)


project(pykingas)
set(TARGET libpykingas)

# Add a library using FindPython's tooling (pybind11 also provides a helper like this)
python_add_library(${TARGET} MODULE ${SOURCES} WITH_SOABI)
target_include_directories(${TARGET} PRIVATE ${pybind11_HEADERS})

# The install directory is the output (wheel) directory
install(TARGETS ${TARGET} DESTINATION ${CMAKE_SOURCE_DIR}/pykingas)