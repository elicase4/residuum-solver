# executables have their own target to be installed, tested, and linked independently

cmake_minimum_required(VERSION 3.18)

find_package(yaml-cpp REQUIRED)

add_library(residuum_heateq_app
	HeatDispatcher.cpp
	parser/BoundaryConditionConfigParser.cpp
	parser/ConductivityConfigParser.cpp
	parser/HeatConfigParser.cpp
	parser/InitialConditionConfigParser.cpp
	parser/MonitorConfigParser.cpp
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

add_executable(heateq
	HeatApplication.cpp
)

target_link_libraries(heateq
    PRIVATE
        residuum_heateq_app
)

# install targets
install(TARGETS heateq
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
