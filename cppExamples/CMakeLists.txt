cmake_minimum_required(VERSION 3.16)
cmake_policy(VERSION 3.12)
project(cppExamples LANGUAGES C CXX)

set(KineticGas_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
find_package(KineticGas REQUIRED)
if (NOT KINETICGAS_INSTALLED)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/.. kineticgas)
endif()

add_executable(basic basic.cpp)
target_link_libraries(basic PRIVATE kineticgas)
message("KINETICGAS_INCLUDE IS : ${KINETICGAS_INCLUDE}")
target_include_directories(basic PRIVATE KINETICGAS_INCLUDE)
