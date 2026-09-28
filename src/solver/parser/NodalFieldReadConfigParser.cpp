#include "solver/parser/NodalFieldReadConfigParser.hpp"

#include <stdexcept>

#include "io/YAMLReader.hpp"

residuum::solver::config::NodalFieldReadConfig::Mode residuum::solver::parser::NodalFieldReadConfigParser::parseMode(const std::string& str) {

	if (str == "expression") {
		return residuum::solver::config::NodalFieldReadConfig::Mode::Expression;
	}

	if (str == "file") {
		return residuum::solver::config::NodalFieldReadConfig::Mode::File;
	}

	throw std::runtime_error("Unknown read mode: " + str);

}

residuum::solver::config::NodalFieldReadConfig residuum::solver::parser::NodalFieldReadConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	residuum::solver::config::NodalFieldReadConfig cfg;

	cfg.mode = NodalFieldReadConfigParser::parseMode(YAMLReader::optional<std::string>(node, "mode", "expression"));

	if (cfg.mode == residuum::solver::config::NodalFieldReadConfig::Mode::File) {
		cfg.file = YAMLReader::required<std::string>(node, "file");
	} else {
		cfg.expression = YAMLReader::required<std::string>(node, "expression");
	}

	cfg.unit = YAMLReader::required<std::string>(node, "unit");

	return cfg;

}
