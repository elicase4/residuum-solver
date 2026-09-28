#include "application/mesh/MeshApplication.hpp"
#include "application/mesh/MeshConfigParser.hpp"
#include "application/mesh/MeshDispatcher.hpp"

#include "solver/logging/LoggerFactory.hpp"
#include "utils/logging/Banner.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {

	residuum::utils::logging::printStartupBanner("mesh", "mesh generation / import");

	if (argc < 2) {
		std::cerr << "Usage: mesh <config.yaml>\n";
		return EXIT_FAILURE;
	}

	// no LoggingConfig exists until the config actually parses, so this is a plain cerr boundary
	residuum::application::mesh::MeshConfig cfg;

	try {
		cfg = residuum::application::mesh::MeshConfigParser::read(argv[1]);
	} catch (const std::exception& e) {
		std::cerr << "[mesh] fatal error: " << e.what() << "\n";
		return EXIT_FAILURE;
	}

	const auto logger = residuum::solver::logging::makeDriverLogger(cfg.logging.driver, "mesh");

	try {

		residuum::application::mesh::MeshApplication app(cfg);
		return app.run();

	} catch (const std::exception& e) {

		logger.error(e.what());
		return EXIT_FAILURE;

	}

}

residuum::application::mesh::MeshApplication::MeshApplication(const MeshConfig& config) : config_(config) {}

int residuum::application::mesh::MeshApplication::run() {
	const bool ok = residuum::application::mesh::MeshDispatcher::run(config_);
	return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
