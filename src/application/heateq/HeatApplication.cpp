#include "application/heateq/HeatApplication.hpp"
#include "application/heateq/HeatDispatcher.hpp"
#include "application/heateq/parser/HeatConfigParser.hpp"

#include "application/cli/AppUsage.hpp"

#include "solver/logging/LoggerFactory.hpp"
#include "utils/logging/Banner.hpp"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

	constexpr const char* kDescription = "heat equation solver";

	// output/log/monitor files are all opened with std::ofstream, which fails silently-then-throws
	// if its parent directory doesn't exist -- create every directory a fresh clone won't have yet
	void ensureOutputDirectories(const residuum::application::heateq::config::HeatConfig& cfg) {

		residuum::solver::logging::ensureLoggingDirectories(cfg.logging);

		if (cfg.output) std::filesystem::create_directories(cfg.output->directory);

		for (const auto& monitor : cfg.monitors) {
			residuum::solver::logging::createParentDirectory(monitor.output.file);
		}

	}

} // namespace

int residuum::application::heateq::runApp(int argc, char** argv) {

	if (argc >= 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
		residuum::application::cli::printAppUsage("heateq", kDescription);
		return EXIT_SUCCESS;
	}

	residuum::utils::logging::printStartupBanner("heateq", kDescription);

	if (argc < 2) {
		residuum::application::cli::printAppUsage("heateq", kDescription);
		return EXIT_FAILURE;
	}

	// no LoggingConfig exists until the config actually parses, so this is a plain cerr boundary
	residuum::application::heateq::config::HeatConfig cfg;

	try {
		cfg = residuum::application::heateq::parser::HeatConfigParser::read(argv[1]);
	} catch (const std::exception& e) {
		std::cerr << "[heateq] fatal error: " << e.what() << "\n";
		return EXIT_FAILURE;
	}

	try {

		ensureOutputDirectories(cfg);

		const auto logger = residuum::solver::logging::makeDriverLogger(cfg.logging.driver, "heateq");

		try {

			residuum::application::heateq::HeatApplication app(cfg);
			return app.run();

		} catch (const std::exception& e) {

			logger.error(e.what());
			return EXIT_FAILURE;

		}

	} catch (const std::exception& e) {

		std::cerr << "[heateq] fatal error: " << e.what() << "\n";
		return EXIT_FAILURE;

	}

}

residuum::application::heateq::HeatApplication::HeatApplication(const residuum::application::heateq::config::HeatConfig& config) : config_(config) {}

int residuum::application::heateq::HeatApplication::run() {

	bool success = false;

	success = residuum::application::heateq::HeatDispatcher::run(config_);

	return success ? EXIT_SUCCESS : EXIT_FAILURE;

}
