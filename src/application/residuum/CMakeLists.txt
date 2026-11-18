cmake_minimum_required(VERSION 3.18)

# The single installed CLI entry point; dispatches to each physics/mesh app's own runApp()
# (see application/cli/CliApp.hpp). Adding a new physics module means adding its library here
# and one more entry in the subcommands vector in ResiduumApplication.cpp.
add_executable(residuum
	ResiduumApplication.cpp
)

target_link_libraries(residuum
	PRIVATE
		residuum_mesh_app
		residuum_heateq_app
)

install(TARGETS residuum
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
