cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
project(wzalgorithm VERSION 0.8.2 LANGUAGES CXX)

set_property(GLOBAL PROPERTY CXX_STANDARD 14)

include(FindSWIG)
include(FindBoost)
include(GNUInstallDirs)
include(CheckIncludeFileCXX)


find_package(GTest)
find_package(SWIG 3.0)
find_package(Boost 1.58.0 REQUIRED COMPONENTS random)
check_include_file_cxx("benchmark/benchmark.h" BENCHMARK_FOUND)


configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/VERSION.in"
    "${CMAKE_CURRENT_BINARY_DIR}/VERSION")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/src/config.h")

include_directories("${CMAKE_CURRENT_BINARY_DIR}/src")


add_custom_target(srctar
    COMMAND tar
        "-C" "${CMAKE_CURRENT_SOURCE_DIR}"
        "-vczf" "libwzalgorithm-${wzalgorithm_VERSION}.tar.gz"
        "--xform='s,^,libwzalgorithm-${wzalgorithm_VERSION}/,'"
        "--exclude-vcs" "--exclude-backups"
        "--exclude=Makefile" "--exclude='*.so*'"
        "src" "test" "CMakeLists.txt"
        "README.md" "COPYING" "VERSION.in" "wzalgorithm.spec.in")


add_subdirectory(src)

if (${GTEST_FOUND})
    include(CTest)
    enable_testing()
    add_subdirectory(test)
else()
    message(AUTHOR_WARNING,
        "GTest was not found; no unit tests will be available and `make check' will not work.")
endif()
