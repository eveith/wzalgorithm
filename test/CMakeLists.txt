add_subdirectory(ann)
add_subdirectory(model)
add_subdirectory(agent)
add_subdirectory(simulation)
file(
    COPY simulation/mock
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/simulation")
