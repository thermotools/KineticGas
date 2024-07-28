# This build assumes that the KineticGas C++ library has already been built and installed with
#   cmake -Dpurecpp=ON ..
#   make install
# Note that the variable KineticGas_DIR is set in this file to help cmake find the installed library.

cmake_minimum_required(VERSION 3.16)
cmake_policy(VERSION 3.12)
project(cppExamples LANGUAGES C CXX)

set(KineticGas_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
find_package(KineticGas REQUIRED) # Assuming that KineticGas has already been installed with `make install`

add_executable(basic basic.cpp)
target_link_libraries(basic PRIVATE kineticgas)

# NOTE: All headers required for KineticGas (including Thermopack, Eigen and json) are automatically added to the include path
#       when we link to the KineticGas library like this.
# If you for some reason want access to the include paths used, they are stored in the environment variable KINETICGAS_INCLUDE
# after you call find_package(KineticGas).

# The path to the dynamic library being linked is in KINETICGAS_LIB.
