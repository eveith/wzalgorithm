cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
PROJECT(winzent-algorithm)

set(winzent-algorithm_MAJOR_VERSION ${WINZENT_MAJOR_VERSION})
set(winzent-algorithm_MINOR_VERSION 1)
set(winzent-algorithm_PATCHLEVEL 0)

set(winzent-algorithm_SOVERSION
    "${winzent-algorithm_MAJOR_VERSION}.${winzent-algorithm_MINOR_VERSION}")

set(winzent-algorithm_SOURCES
    REvol.cpp)

set(winzent-algorithm_HEADERS
    algorithm_global.h
    REvol.h)


# Tweak for C++14:

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU"
        OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    add_compile_options(-std=c++1y)
endif()


add_library(winzent-algorithm SHARED ${winzent-algorithm_SOURCES})
# target_link_libraries(winzent-algorithm ...)
set_target_properties(winzent-algorithm 
    PROPERTIES SOVERSION ${winzent-algorithm_SOVERSION})

install(TARGETS winzent-algorithm LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${winzent-algorithm_HEADERS} DESTINATION include/winzent/algorithm)
