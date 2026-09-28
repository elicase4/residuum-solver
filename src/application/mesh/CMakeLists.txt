cmake_minimum_required(VERSION 3.18)

find_package(yaml-cpp REQUIRED)

# Mesh application library
add_library(residuum_mesh_app
    MeshConfigParser.cpp
    MeshDispatcher.cpp
    GmshMeshGenerator.cpp
)

target_link_libraries(residuum_mesh_app
    PUBLIC
        residuum_headers
        residuum_mesh
        residuum_io
        residuum_solver
        yaml-cpp
)

# Mesh executable
add_executable(mesh
    MeshApplication.cpp
)

target_link_libraries(mesh
    PRIVATE
        residuum_mesh_app
)

# Install target
install(TARGETS mesh
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
