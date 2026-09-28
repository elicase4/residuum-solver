add_library(residuum_mesh
	Mesh.cpp
	exchange/gmsh/IntermediateMesh.cpp
	exchange/gmsh/Utils.cpp
	exchange/gmsh/MeshConverter.cpp
	generator/BlockMesh2D.cpp
	generator/BlockMesh3D.cpp
)

target_link_libraries(residuum_mesh
	PUBLIC
		residuum_headers
)
