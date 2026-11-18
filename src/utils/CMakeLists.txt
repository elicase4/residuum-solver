find_path(EXPRTK_INCLUDE_DIR
    NAMES exprtk.hpp
    PATH_SUFFIXES include
)

if (NOT EXPRTK_INCLUDE_DIR)
    message(FATAL_ERROR "Could not find exprtk.hpp")
endif()

add_library(residuum_expression
	expression/ScalarExpression.cpp
	expression/VectorExpression.cpp
	expression/TensorExpression.cpp
)

target_link_libraries(residuum_expression
	PUBLIC
		residuum_headers
)

# exprtk is header-only and located via find_path rather than find_package, so this absolute path
# is only valid for in-tree consumers; cmake/ResiduumConfig.cmake.in re-runs find_path and adds
# the result to Residuum::residuum_expression's INTERFACE_INCLUDE_DIRECTORIES for installed consumers.
target_include_directories(residuum_expression
    PUBLIC
        $<BUILD_INTERFACE:${EXPRTK_INCLUDE_DIR}>
)
