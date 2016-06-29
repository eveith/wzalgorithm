cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
PROJECT(winzent-algorithm)

set(winzent-algorithm_VERSION_MAJOR ${WINZENT_VERSION_MAJOR})
set(winzent-algorithm_VERSION_MINOR 1)
set(winzent-algorithm_PATCHLEVEL 0)

set(winzent-algorithm_SOVERSION
    "${winzent-algorithm_VERSION_MAJOR}.${winzent-algorithm_VERSION_MINOR}")

set(winzent-algorithm_SOURCES
    REvol.cpp
    SimulatedAnnealing.cpp
    ParticleSwarmOptimization.cpp)

set(winzent-algorithm_HEADERS
    algorithm_global.h
    REvol.h
    SimulatedAnnealing.h
    ParticleSwarmOptimization.h)

set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Qt5Core REQUIRED)


if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU"
        OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funroll-loops")
endif()


add_library(winzent-algorithm SHARED ${winzent-algorithm_SOURCES})
target_include_directories(winzent-algorithm
    PUBLIC ${winzent-common_SOURCE_DIR})
target_link_libraries(winzent-algorithm
    m
    Qt5::Core)
set_target_properties(winzent-algorithm
    PROPERTIES SOVERSION ${winzent-algorithm_SOVERSION}
    CXX_STANDARD 14)

install(TARGETS winzent-algorithm LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${winzent-algorithm_HEADERS} DESTINATION include/winzent/algorithm)
