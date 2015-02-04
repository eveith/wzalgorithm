project(Winzent)
cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

set(WINZENT_VERSION_MAJOR 0)


include(CTest)
include(FindPkgConfig)


add_subdirectory(src)
add_subdirectory(test)
