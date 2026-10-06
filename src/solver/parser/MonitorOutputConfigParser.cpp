#include "solver/parser/MonitorOutputConfigParser.hpp"

#include "io/YAMLReader.hpp"

residuum::solver::config::MonitorOutputConfig residuum::solver::parser::MonitorOutputConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	residuum::solver::config::MonitorOutputConfig cfg;

	cfg.console = YAMLReader::optional<bool>(node, "console", true);
	cfg.file = YAMLReader::optional<std::string>(node, "file", "");

	return cfg;

}
