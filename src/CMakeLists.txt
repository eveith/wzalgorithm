set(wzalgorithm_SOURCES
    REvol.cpp
    Result.cpp
    ParticleSwarmOptimization.cpp)

set(wzalgorithm_HEADERS
    REvol.h
    Result.h
    ParticleSwarmOptimization.h)


set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU"
        OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funroll-loops")
endif()


add_library(wzalgorithm ${wzalgorithm_SOURCES})

target_include_directories(wzalgorithm
    PUBLIC ${Boost_INCLUDE_DIRS})

target_link_libraries(wzalgorithm
    m
    stdc++
    pthread
    ${Boost_LIBRARIES})

set_target_properties(wzalgorithm
    PROPERTIES SOVERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")


SET(PKG_CONFIG_REQUIRES "Boost")
SET(PKG_CONFIG_LIBDIR "\${prefix}/${CMAKE_INSTALL_LIBDIR}")
SET(PKG_CONFIG_INCLUDEDIR "\${prefix}/include")
SET(PKG_CONFIG_LIBS "-L\${libdir} ${wzalgorithm_LIBRARIES} -lwzalgorithm")
SET(PKG_CONFIG_CFLAGS "-I\${includedir}")

CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/libwzalgorithm.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/libwzalgorithm.pc)


#
# Bindings for other Languages
#


if (SWIG_FOUND)
    if (PYTHONLIBS_FOUND)
        set_source_files_properties(wzalgorithm-python.i PROPERTIES
            CPLUSPLUS ON
            SWIGFLAGS "-includeall")
        include_directories(${PYTHON_INCLUDE_DIRS})

        if (CMAKE_VERSION VERSION_LESS "3.8.0")
            swig_add_module(wzalgorithm python wzalgorithm-python.i)
        else()
            swig_add_library(wzalgorithm
                LANGUAGE python
                SOURCES wzalgorithm-python.i)
        endif()

        swig_link_libraries(wzalgorithm wzalgorithm ${PYTHON_LIBRARIES})

        install(
            FILES ${CMAKE_CURRENT_BINARY_DIR}/wzalgorithm.py
            DESTINATION ${PYTHON_INSTALL_PACKAGES_DIR})
        install(TARGETS ${SWIG_MODULE_wzalgorithm_REAL_NAME}
            LIBRARY DESTINATION ${PYTHON_INSTALL_PACKAGES_DIR})
    endif()
endif()


#
# Installation
#

install(TARGETS wzalgorithm
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
    FILES
        ${wzalgorithm_HEADERS}
        ${CMAKE_CURRENT_BINARY_DIR}/config.h
    DESTINATION include/wzalgorithm)
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/libwzalgorithm.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
