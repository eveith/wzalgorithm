set(wzalgorithm_SOURCES
    REvol.cpp
    ParticleSwarmOptimization.cpp)

set(wzalgorithm_HEADERS
    config.h
    REvol.h
    ParticleSwarmOptimization.h)


set(CMAKE_INCLUDE_CURRENT_DIR ON)

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
    ${Boost_LIBRARIES})

set_target_properties(wzalgorithm
    PROPERTIES SOVERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
    CXX_STANDARD 14)


install(TARGETS wzalgorithm
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${wzalgorithm_HEADERS} DESTINATION include/winzent/algorithm)
