project(Winzent)
cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

set(WINZENT_VERSION_MAJOR 0)
set(WINZENT_VERSION_MINOR 6)
set(WINZENT_VERSION_PATCH 0)
set(WINZENT_VERSION
    "${WINZENT_VERSION_MAJOR}.${WINZENT_VERSION_MINOR}.${WINZENT_VERSION_PATCH}")


include(CTest)
include(FindPkgConfig)


configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/winzent.spec.in"
    "${CMAKE_CURRENT_BINARY_DIR}/winzent.spec")

add_custom_target(srctar
    ALL
    COMMAND tar 
        "-C" "${CMAKE_CURRENT_SOURCE_DIR}"
        "-vcJf" "winzent-${WINZENT_VERSION}.tar.xz"
        "--exclude-vcs" "--exclude-backups"
        "--exclude=Makefile" "--exclude=src/flavormap" "--exclude='*.so*'"
        "src" "test" "CMakeLists.txt" "README")


add_subdirectory(src)
add_subdirectory(test)
