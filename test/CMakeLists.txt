set(test-wzalgorithm_SOURCES
    REvolTest.cpp
    ParticleSwarmOptimizationTest.cpp
    test_wzalgorithm.cpp)

set(test-wzalgorithm_HEADERS
    REvolTest.h
    ParticleSwarmOptimizationTest.h)


set(benchmark-wzalgorithm_SOURCES
    REvolBenchmark.cpp
    benchmark_wzalgorithm.cpp)

set(benchmark-wzalgorithm_HEADERS
    REvolBenchmark.h)


set(CMAKE_INCLUDE_CURRENT_DIR ON)


add_executable(test_wzalgorithm
    ${test-wzalgorithm_SOURCES})

add_test(test_wzalgorithm
    test_wzalgorithm)

target_include_directories(test_wzalgorithm
    PUBLIC "${wzalgorithm_SOURCE_DIR}/src"
    PUBLIC ${Boost_INCLUDE_DIRS}
    PUBLIC ${GTEST_INCLUDE_DIRS})

target_link_libraries(test_wzalgorithm
    m
    stdc++
    ${Boost_LIBRARIES}
    ${GTEST_LIBRARIES}
    wzalgorithm)

set_target_properties(test_wzalgorithm
    PROPERTIES CXX_STANDARD 14)


if (${BENCHMARK_FOUND})
    add_executable(benchmark_wzalgorithm
        ${benchmark-wzalgorithm_SOURCES})

    add_test(benchmark_wzalgorithm
        benchmark_wzalgorithm)

    target_include_directories(benchmark_wzalgorithm
        PUBLIC "${wzalgorithm_SOURCE_DIR}/src"
        PUBLIC ${Boost_INCLUDE_DIRS}
        PUBLIC ${GTEST_INCLUDE_DIRS})

    target_link_libraries(benchmark_wzalgorithm
        m
        stdc++
        ${Boost_LIBRARIES}
        benchmark
        wzalgorithm)

    set_target_properties(benchmark_wzalgorithm
        PROPERTIES CXX_STANDARD 14)
endif()
