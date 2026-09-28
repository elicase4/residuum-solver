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
# won't re-resolve when the exported target is consumed on a different machine.
target_include_directories(residuum_expression
    PUBLIC
        ${EXPRTK_INCLUDE_DIR}
)
