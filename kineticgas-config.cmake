# Checks whether KineticGas has been built and installed.
# Result variables:
#   KINETICGAS_FOUND : TRUE
#   KINETICGAS_INSTALLED : Whether the dynamic library is found where expected (TRUE/FALSE).
#   KINETICGAS_LIB : Path to the dynamic library
#   KINETICGAS_INCLUDE : List of include paths to KineticGas headers and all dependencies
#   KINETICGAS_INCLUDE_ONLY : Include path to KineticGas headers


set(KINETICGAS_LIB ${CMAKE_CURRENT_LIST_DIR}/lib/libkineticgas${CMAKE_SHARED_LIBRARY_SUFFIX})

if(NOT EXISTS ${KINETICGAS_LIB})
    unset(KINETICGAS_LIB)
    set(KINETICGAS_INSTALLED FALSE)
    return()
endif()

set(KINETICGAS_FOUND TRUE)
set(KINETICGAS_INSTALLED TRUE)
set(KINETICGAS_INCLUDE_ONLY ${CMAKE_CURRENT_LIST_DIR}/cpp)

if(NOT TARGET kineticgas)

    find_package(thermopack REQUIRED)
    if(NOT THERMOPACK_INSTALLED)
        add_subdirectory(${THERMOPACK_ROOT}) # We need to build thermopack, because it has not been installed.
    endif()

    set(JSON_INCLUDE ${CMAKE_CURRENT_LIST_DIR}/cpp/external/json/include/) # The directory json/include/nlohmann/ is here (https://github.com/nlohmann/json)
    set(EIGEN_INCLUDE ${CMAKE_CURRENT_LIST_DIR}/cpp/external/eigen/)

    add_library(kineticgas SHARED IMPORTED)

    set(ALL_KINETICGAS_INCLUDE_DIRS ${KINETICGAS_INCLUDE_ONLY} ${THERMOPACK_INCLUDE} ${JSON_INCLUDE} ${EIGEN_INCLUDE})
    list(JOIN ALL_KINETICGAS_INCLUDE_DIRS ";" KINETICGAS_INCLUDE)

    set_target_properties(kineticgas PROPERTIES 
                            IMPORTED_LOCATION ${KINETICGAS_LIB}
                            INTERFACE_INCLUDE_DIRECTORIES "${KINETICGAS_INCLUDE}")
    target_link_libraries(kineticgas INTERFACE thermopack)
    target_compile_definitions(kineticgas INTERFACE NOPYTHON)

    


endif()