#include "application/cli/CliApp.hpp"
#include "application/cli/Subcommand.hpp"

#include "application/heateq/HeatApplication.hpp"
#include "application/mesh/MeshApplication.hpp"

#include <vector>

int main(int argc, char** argv) {

	const std::vector<residuum::application::cli::Subcommand> subcommands = {
		{"mesh", "mesh generation / import", &residuum::application::mesh::runApp},
		{"heateq", "heat equation solver", &residuum::application::heateq::runApp},
	};

	return residuum::application::cli::runSubcommandDispatch(argc, argv, subcommands);

}
