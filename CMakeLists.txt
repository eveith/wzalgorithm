cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
project(wzalgorithm VERSION 0.8.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

include(FindSWIG)
include(FindBoost)
include(GNUInstallDirs)
include(CheckIncludeFileCXX)


find_package(GTest)
find_package(SWIG 3.0)
find_package(Boost 1.58.0 REQUIRED COMPONENTS random)
check_include_file_cxx("benchmark/benchmark.h" BENCHMARK_FOUND)

if (SWIG_FOUND)
    include(UseSWIG)

    include(FindPythonLibs)
    find_package(PythonLibs)
endif()


if (PYTHONLIBS_FOUND)
    execute_process(
        COMMAND python3
            -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
        OUTPUT_VARIABLE PYTHON_INSTALL_PACKAGES_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

set(PYTHON_INSTALL_PACKAGES_DIR ${CMAKE_INSTALL_LIBDIR}
    CACHE PATH "Installation path for the Python bindings")


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
        "--exclude=__pycache__" "--exclude='*.pyc'"
        "--exclude='CMakeCache*'" "--exclude=CMakeFiles"
        "src" "test" "CMakeLists.txt"
        "README.md" "COPYING" "VERSION.in" "wzalgorithm.spec.in")


add_subdirectory(src)
include(CTest)
enable_testing()
add_subdirectory(test)
