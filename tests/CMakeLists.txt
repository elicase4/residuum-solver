# Add unit tests
function(add_residuum_test TEST_NAME SOURCE_FILE)
	add_executable(${TEST_NAME} ${SOURCE_FILE})

	# Data directories
	set(TEST_DATA_DIR "${CMAKE_SOURCE_DIR}/tests/data")
	set(TEST_OUTPUT_DIR "${CMAKE_BINARY_DIR}/test_outputs/${TEST_NAME}")
	file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

	# Organize binaries
	set_target_properties(${TEST_NAME} PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/tests"
	)
	
	# Cpp macro for test paths
	target_compile_definitions(${TEST_NAME} PRIVATE
		TEST_OUTPUT_PATH="${TEST_OUTPUT_DIR}"
		TEST_DATA_PATH="${TEST_DATA_DIR}"
	)

	# Link against GTest libraries
	target_link_libraries(${TEST_NAME}
		PRIVATE
			residuum_headers
			${ARGN}
			GTest::gtest_main
	)

	# Discover unit tests automatically
	gtest_discover_tests(${TEST_NAME}
		PROPERTIES TIMEOUT 30
	)

endfunction()

add_residuum_test(UnitTest_Mesh_Generator_BlockMesh2D
	unit/mesh/generator/BlockMesh2D.cpp
	residuum_mesh residuum_io
)

add_residuum_test(UnitTest_Mesh_Generator_BlockMesh3D
	unit/mesh/generator/BlockMesh3D.cpp
	residuum_mesh residuum_io
)

add_residuum_test(UnitTest_Mesh_Exchange_Gmsh_MeshConverter
	unit/mesh/exchange/gmsh/MeshConverter.cpp
	residuum_mesh residuum_io
)

add_residuum_test(UnitTest_Fem_Basis_Lagrange
	unit/fem/basis/Lagrange.cpp
)

add_residuum_test(UnitTest_Fem_Geometry_JacobianTransform
	unit/fem/geometry/JacobianTransform.cpp
)

add_residuum_test(UnitTest_Fem_Quadrature_GaussQuadrature
	unit/fem/quadrature/GaussQuadrature.cpp
)

add_residuum_test(UnitTest_Fem_Dispatch_DiscretizationDispatch
	unit/fem/dispatch/DiscretizationDispatch.cpp
)

add_residuum_test(UnitTest_Fem_Quantity_QuantityUnits
	unit/fem/quantity/QuantityUnits.cpp
)

add_residuum_test(UnitTest_Fem_Quantity_MonitorGroup
	unit/fem/quantity/MonitorGroup.cpp
	residuum_mesh residuum_io
)

add_residuum_test(UnitTest_Linalg_Operations_VectorOps
	unit/linalg/operations/VectorOps.cpp
)

add_residuum_test(UnitTest_Linalg_Operations_MatrixOps
	unit/linalg/operations/MatrixOps.cpp
)

add_residuum_test(UnitTest_Linalg_Solver_Iterative_Cg_Solver
	unit/linalg/solver/iterative/cg/Solver.cpp
)

add_residuum_test(UnitTest_Linalg_Solver_Iterative_Cg_CGRunner
	unit/linalg/solver/iterative/cg/CGRunner.cpp
)

add_residuum_test(UnitTest_Linalg_Solver_Iterative_Gmres_Solver
	unit/linalg/solver/iterative/gmres/Solver.cpp
)

add_residuum_test(UnitTest_Linalg_Solver_Iterative_Gmres_GMRESRunner
	unit/linalg/solver/iterative/gmres/GMRESRunner.cpp
)

add_residuum_test(UnitTest_Linalg_Solver_Preconditioner_Jacobi
	unit/linalg/solver/preconditioner/Jacobi.cpp
)

add_residuum_test(UnitTest_Solver_Linear_LinearSolverFactory
	unit/solver/linear/LinearSolverFactory.cpp
)

add_residuum_test(UnitTest_Solver_Nonlinear_Newton
	unit/solver/nonlinear/Newton.cpp
)

add_residuum_test(UnitTest_IO_Fieldio_FieldIO
	unit/io/fieldio/FieldIO.cpp
	residuum_mesh residuum_io
)

add_residuum_test(UnitTest_IO_Visualization_VisualizationWriter
	unit/io/visualization/VisualizationWriter.cpp
	residuum_mesh residuum_io
)

add_residuum_test(UnitTest_IO_Fieldio_TimeSeriesValueSource
	unit/io/fieldio/TimeSeriesValueSource.cpp
	residuum_mesh residuum_io
)

add_residuum_test(UnitTest_IO_MeshIO
	unit/io/MeshIO.cpp
	residuum_mesh residuum_io
)

add_residuum_test(UnitTest_IO_GmshReader
	unit/io/GmshReader.cpp
	residuum_mesh residuum_io
)

add_residuum_test(UnitTest_Fem_Form_FormRegistry
	unit/fem/form/FormRegistry.cpp
	residuum_expression
)

add_residuum_test(UnitTest_Equation_Heateq_ScalarExpressionWrappers
	unit/equation/heateq/ScalarExpressionWrappers.cpp
	residuum_expression
)

add_residuum_test(UnitTest_Equation_Heateq_NodalDataFormEquivalence
	unit/equation/heateq/NodalDataFormEquivalence.cpp
	residuum_mesh residuum_io
)

add_residuum_test(UnitTest_Equation_Heateq_QuantityEvaluator
	unit/equation/heateq/QuantityEvaluator.cpp
	residuum_mesh residuum_io
)

add_residuum_test(UnitTest_Equation_Heateq_Evaluator_ConductivityModel
	unit/equation/heateq/evaluator/ConductivityModel.cpp
	residuum_expression
)

add_residuum_test(UnitTest_Equation_Heateq_Evaluator_SpecificHeatModel
	unit/equation/heateq/evaluator/SpecificHeatModel.cpp
	residuum_expression
)

add_residuum_test(UnitTest_Application_Heateq_Problem_HeatProblem
	unit/application/heateq/problem/HeatProblem.cpp
	residuum_mesh residuum_io residuum_solver residuum_expression
)

add_residuum_test(UnitTest_Solver_Stage_SteadyStage
	unit/solver/stage/SteadyStage.cpp
	residuum_mesh residuum_io residuum_solver residuum_expression
)

add_residuum_test(UnitTest_Solver_Stage_BackwardEulerStage
	unit/solver/stage/BackwardEulerStage.cpp
	residuum_mesh residuum_io residuum_solver residuum_expression
)

add_residuum_test(IntegrationTest_Fem_Heateq_Cpu_HeatEquationMinimal
	integration/fem/heateq/cpu/HeatEquationMinimal2D.cpp
	residuum_mesh
)

add_residuum_test(IntegrationTest_Fem_Heateq_Cpu_HeatEquationMinimal3D
	integration/fem/heateq/cpu/HeatEquationMinimal3D.cpp
	residuum_mesh
)

add_residuum_test(IntegrationTest_Full_Heateq_Cpu_HeatEquationMinimal
	integration/full/heateq/cpu/HeatEquationMinimal2D.cpp
	residuum_mesh residuum_io
)

add_residuum_test(IntegrationTest_Full_Heateq_Cpu_HeatEquationMinimal3D
	integration/full/heateq/cpu/HeatEquationMinimal3D.cpp
	residuum_mesh residuum_io
)

add_residuum_test(IntegrationTest_Full_Heateq_Cpu_HeatEquationFileMode
	integration/full/heateq/cpu/HeatEquationFileMode2D.cpp
	residuum_mesh residuum_io residuum_solver residuum_heateq_app
)

add_residuum_test(IntegrationTest_Full_Heateq_Cpu_HeatEquationTransientBackwardEuler2D
	integration/full/heateq/cpu/HeatEquationTransientBackwardEuler2D.cpp
	residuum_mesh residuum_io residuum_solver residuum_expression
)

add_residuum_test(IntegrationTest_Full_Heateq_Cpu_HeatEquationSteadyNonlinearConductivity2D
	integration/full/heateq/cpu/HeatEquationSteadyNonlinearConductivity2D.cpp
	residuum_mesh residuum_io residuum_solver residuum_expression
)
