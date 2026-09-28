find_package(yaml-cpp REQUIRED)

add_library(residuum_io
	MeshIO.cpp
	GmshReader.cpp
	visualization/VTKWriter.cpp
	visualization/VTUWriter.cpp
	visualization/VTUSeriesWriter.cpp
	visualization/VisualizationWriter.cpp
	YAMLReader.cpp
)

target_link_libraries(residuum_io
	PUBLIC
		residuum_headers
		residuum_mesh
		yaml-cpp
)
