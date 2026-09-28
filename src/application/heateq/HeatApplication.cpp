#include "application/heateq/HeatApplication.hpp"
#include "application/heateq/HeatDispatcher.hpp"
#include "application/heateq/parser/HeatConfigParser.hpp"

#include "solver/logging/LoggerFactory.hpp"
#include "utils/logging/Banner.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {

	residuum::utils::logging::printStartupBanner("heateq", "heat equation solver");

	if (argc < 2) {
		std::cerr << "Usage: heateq <config.yaml>\n";
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

	const auto logger = residuum::solver::logging::makeDriverLogger(cfg.logging.driver, "heateq");

	try {

		residuum::application::heateq::HeatApplication app(cfg);
		return app.run();

	} catch (const std::exception& e) {

		logger.error(e.what());
		return EXIT_FAILURE;

	}

}

residuum::application::heateq::HeatApplication::HeatApplication(const residuum::application::heateq::config::HeatConfig& config) : config_(config) {}

int residuum::application::heateq::HeatApplication::run() {

	bool success = false;

	success = residuum::application::heateq::HeatDispatcher::run(config_);

	return success ? EXIT_SUCCESS : EXIT_FAILURE;

}
