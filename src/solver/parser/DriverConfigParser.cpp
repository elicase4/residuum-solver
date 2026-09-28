#include "solver/parser/DriverConfigParser.hpp"
#include "io/YAMLReader.hpp"

residuum::solver::config::DriverConfig::Type residuum::solver::parser::DriverConfigParser::parseDriverType(const std::string& str) {

	if (str == "steady")
		return config::DriverConfig::Type::Steady;
	if (str == "transient")
		return config::DriverConfig::Type::Transient;
	if (str == "pseudo_transient")
		return config::DriverConfig::Type::PseudoTransient;

	throw std::runtime_error("Unknown driver type: '" + str + "'. Valid options: steady, transient, pseudo_transient");

}

residuum::solver::config::DriverConfig residuum::solver::parser::DriverConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	config::DriverConfig cfg;

	cfg.type = residuum::solver::parser::DriverConfigParser::parseDriverType(YAMLReader::required<std::string>(node, "type"));

	return cfg;

}
