project(test-winzent-algorithm CXX)
cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)


set(test-winzent-algorithm_SOURCES
    REvolTest.cpp
    tst_algorithm.cpp)

set(test-winzent-algorithm_HEADERS
    REvolTest.h)


set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Qt5Core REQUIRED)
find_package(Qt5Test REQUIRED)
find_package(Boost "1.56.0" REQUIRED COMPONENTS random )
pkg_check_modules(LOG4CXX REQUIRED liblog4cxx>=0.10.0)


# Tweak for C++11:

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU"
        OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    add_compile_options(-std=c++1y)
endif()



add_executable(tst_algorithm
    ${test-winzent-algorithm_SOURCES})
add_test(tst_algorithm
    tst_algorithm)
target_include_directories(tst_algorithm
    PUBLIC ${winzent-algorithm_SOURCE_DIR}
    PUBLIC ${winzent-test-common_SOURCE_DIR}
    PUBLIC ${BOOST_INCLUDE_DIRS}
    PUBLIC ${LOG4CXX_INCLUDE_DIRS})
target_link_libraries(tst_algorithm
    m
    stdc++
    Qt5::Core
    Qt5::Test
    ${LOG4CXX_LIBRARIES}
    ${BOOST_LIBRARIES}
    winzent-test-common
    winzent-algorithm)

