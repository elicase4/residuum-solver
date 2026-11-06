#include "application/mesh/MeshApplication.hpp"
#include "application/mesh/MeshConfigParser.hpp"
#include "application/mesh/MeshDispatcher.hpp"

#include "application/cli/AppUsage.hpp"

#include "solver/logging/LoggerFactory.hpp"
#include "utils/logging/Banner.hpp"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

	constexpr const char* kDescription = "mesh generation / import";

	// output/log files are all opened with std::ofstream, which fails silently-then-throws if its
	// parent directory doesn't exist -- create every directory a fresh clone won't have yet
	void ensureOutputDirectories(const residuum::application::mesh::MeshConfig& cfg) {

		residuum::solver::logging::ensureLoggingDirectories(cfg.logging);
		residuum::solver::logging::createParentDirectory(cfg.outputFile);

	}

} // namespace

int residuum::application::mesh::runApp(int argc, char** argv) {

	if (argc >= 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
		residuum::application::cli::printAppUsage("mesh", kDescription);
		return EXIT_SUCCESS;
	}

	residuum::utils::logging::printStartupBanner("mesh", kDescription);

	if (argc < 2) {
		residuum::application::cli::printAppUsage("mesh", kDescription);
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

	try {

		ensureOutputDirectories(cfg);

		const auto logger = residuum::solver::logging::makeDriverLogger(cfg.logging.driver, "mesh");

		try {

			residuum::application::mesh::MeshApplication app(cfg);
			return app.run();

		} catch (const std::exception& e) {

			logger.error(e.what());
			return EXIT_FAILURE;

		}

	} catch (const std::exception& e) {

		std::cerr << "[mesh] fatal error: " << e.what() << "\n";
		return EXIT_FAILURE;

	}

}

residuum::application::mesh::MeshApplication::MeshApplication(const MeshConfig& config) : config_(config) {}

int residuum::application::mesh::MeshApplication::run() {
	const bool ok = residuum::application::mesh::MeshDispatcher::run(config_);
	return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
