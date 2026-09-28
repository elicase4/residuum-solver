find_package(yaml-cpp REQUIRED)

add_library(residuum_solver
	SolverInstance.cpp
	parser/BackendConfigParser.cpp
	parser/DriverConfigParser.cpp
	parser/DiscretizationConfigParser.cpp
	parser/MeshConfigParser.cpp
	parser/SolverConfigParser.cpp
	parser/LinearSolverConfigParser.cpp
	parser/NonlinearSolverConfigParser.cpp
	parser/TimeStepperConfigParser.cpp
	parser/TimeStepSizeConfigParser.cpp
	parser/OutputConfigParser.cpp
	parser/LoggingConfigParser.cpp
	parser/NodalFieldReadConfigParser.cpp
)

target_link_libraries(residuum_solver
	PUBLIC
		residuum_headers
		residuum_expression
		yaml-cpp
)
