# Add benchmark executables
function(add_residuum_benchmark BENCHMARK_NAME SOURCE_FILE)
	add_executable(${BENCHMARK_NAME} ${SOURCE_FILE})

	# Results directory
	set(BENCHMARK_RESULTS_DIR "${CMAKE_BINARY_DIR}/benchmark_results")
	file(MAKE_DIRECTORY "${BENCHMARK_RESULTS_DIR}")

	# Organize binaries
	set_target_properties(${BENCHMARK_NAME} PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/benchmarks"
	)

	# Cpp macro for the CSV output path
	target_compile_definitions(${BENCHMARK_NAME} PRIVATE
		BENCHMARK_RESULTS_PATH="${BENCHMARK_RESULTS_DIR}"
	)

	target_link_libraries(${BENCHMARK_NAME}
		PRIVATE
			residuum_headers
			${ARGN}
	)

	# not registered with ctest -- these are measurements, not pass/fail checks

endfunction()

add_subdirectory(linalg)
add_subdirectory(equation)
