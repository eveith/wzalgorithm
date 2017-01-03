cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
project(wzalgorithm VERSION 0.7.0 LANGUAGES CXX)


include(GNUInstallDirs)


find_package(GTest)
find_package(Boost 1.58.0 REQUIRED COMPONENTS random)


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
        "-vcJf" "wzalgorithm-${WZALGORITHM_VERSION}.tar.xz"
        "--xform='s,^,wzalgorithm-${WZALGORITHM_VERSION}/,'"
        "--exclude-vcs" "--exclude-backups"
        "--exclude=Makefile" "--exclude='*.so*'"
        "src" "test" "CMakeLists.txt"
        "README" "COPYING" "VERSION.in" "wzalgorithm.spec.in")


add_subdirectory(src)

if (${GTEST_FOUND})
    include(CTest)
    enable_testing()
    add_subdirectory(test)
else()
    message(AUTHOR_WARNING,
        "GTest was not found; no unit tests will be available and `make check' will not work.")
endif()
