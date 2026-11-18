cmake_minimum_required(VERSION 3.18)

find_package(yaml-cpp REQUIRED)

# Mesh application library; MeshApplication.cpp's runApp() is the entry point the residuum
# subcommand dispatcher calls -- there is no standalone mesh executable (see application/residuum)
add_library(residuum_mesh_app
    MeshApplication.cpp
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
