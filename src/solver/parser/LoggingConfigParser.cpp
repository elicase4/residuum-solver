#include "solver/parser/LoggingConfigParser.hpp"

#include <stdexcept>

#include "io/YAMLReader.hpp"

residuum::solver::config::LoggerConfig::Type residuum::solver::parser::LoggingConfigParser::parseLoggerType(const std::string& str) {

	if (str == "console") {
		return residuum::solver::config::LoggerConfig::Type::Console;
	}

	if (str == "none") {
		return residuum::solver::config::LoggerConfig::Type::None;
	}

	throw std::runtime_error("Unknown logger type: " + str);

}

residuum::solver::config::LoggingConfig residuum::solver::parser::LoggingConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	residuum::solver::config::LoggingConfig cfg;

	if (!node) {
		return cfg;
	}

	const YAML::Node& linearNode = node["linear"];
	if (linearNode) {
		cfg.linear.type = LoggingConfigParser::parseLoggerType(YAMLReader::optional<std::string>(linearNode, "type", "console"));
		cfg.linear.textFile = YAMLReader::optional<std::string>(linearNode, "text_file", "");
		cfg.linear.csvFile = YAMLReader::optional<std::string>(linearNode, "csv_file", "");
		cfg.linear.interval = YAMLReader::optional<Index>(linearNode, "interval", 1);
	}

	const YAML::Node& timestepperNode = node["timestepper"];
	if (timestepperNode) {
		cfg.timestepper.type = LoggingConfigParser::parseLoggerType(YAMLReader::optional<std::string>(timestepperNode, "type", "console"));
		cfg.timestepper.textFile = YAMLReader::optional<std::string>(timestepperNode, "text_file", "");
		cfg.timestepper.csvFile = YAMLReader::optional<std::string>(timestepperNode, "csv_file", "");
		cfg.timestepper.interval = YAMLReader::optional<Index>(timestepperNode, "interval", 1);
	}

	const YAML::Node& nonlinearNode = node["nonlinear"];
	if (nonlinearNode) {
		cfg.nonlinear.type = LoggingConfigParser::parseLoggerType(YAMLReader::optional<std::string>(nonlinearNode, "type", "console"));
		cfg.nonlinear.textFile = YAMLReader::optional<std::string>(nonlinearNode, "text_file", "");
		cfg.nonlinear.csvFile = YAMLReader::optional<std::string>(nonlinearNode, "csv_file", "");
		cfg.nonlinear.interval = YAMLReader::optional<Index>(nonlinearNode, "interval", 1);
	}

	const YAML::Node& driverNode = node["driver"];
	if (driverNode) {
		cfg.driver.type = LoggingConfigParser::parseLoggerType(YAMLReader::optional<std::string>(driverNode, "type", "console"));
		cfg.driver.textFile = YAMLReader::optional<std::string>(driverNode, "text_file", "");
	}

	return cfg;

}
