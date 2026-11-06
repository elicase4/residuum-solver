cmake_minimum_required(VERSION 3.18)

find_package(yaml-cpp REQUIRED)

# HeatApplication.cpp's runApp() is the entry point the residuum subcommand dispatcher calls --
# there is no standalone heateq executable (see application/residuum)
add_library(residuum_heateq_app
	HeatApplication.cpp
	HeatDispatcher.cpp
	parser/BoundaryConditionConfigParser.cpp
	parser/ConductivityConfigParser.cpp
	parser/HeatConfigParser.cpp
	parser/InitialConditionConfigParser.cpp
	parser/MonitorConfigParser.cpp
	parser/OutputFieldConfigParser.cpp
	parser/ScalarMaterialPropertyConfigParser.cpp
	parser/SourceConfigParser.cpp
	parser/SpecificHeatConfigParser.cpp
)

target_link_libraries(residuum_heateq_app
	PUBLIC
		residuum_headers
	PRIVATE
		residuum_expression
		residuum_io
		residuum_solver
		yaml-cpp
)
