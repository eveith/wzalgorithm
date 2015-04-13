project(Winzent)
cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

set(WINZENT_VERSION_MAJOR 0)
set(WINZENT_VERSION_MINOR 6)
set(WINZENT_VERSION_PATCH 0)
set(WINZENT_VERSION
    "${WINZENT_VERSION_MAJOR}.${WINZENT_VERSION_MINOR}.${WINZENT_VERSION_PATCH}")


include(CTest)
include(FindPkgConfig)
include(GNUInstallDirs)


configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/winzent.spec.in"
    "${CMAKE_CURRENT_BINARY_DIR}/winzent.spec")

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/VERSION.in"
    "${CMAKE_CURRENT_BINARY_DIR}/VERSION")


add_custom_target(srctar
    COMMAND tar
        "-C" "${CMAKE_CURRENT_SOURCE_DIR}"
        "-vcJf" "winzent-${WINZENT_VERSION}.tar.xz"
        "--xform='s,^,winzent-${WINZENT_VERSION}/,'"
        "--exclude-vcs" "--exclude-backups"
        "--exclude=Makefile" "--exclude=src/flavormap" "--exclude='*.so*'"
        "src" "test" "CMakeLists.txt"
        "README" "COPYING" "VERSION.in" "winzent.spec.in")


add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(experiment)
