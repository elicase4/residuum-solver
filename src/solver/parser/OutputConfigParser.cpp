#include "solver/parser/OutputConfigParser.hpp"

#include <stdexcept>

#include "io/YAMLReader.hpp"

residuum::solver::config::OutputConfig residuum::solver::parser::OutputConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	residuum::solver::config::OutputConfig cfg;

	cfg.directory = YAMLReader::optional<std::string>(node, "directory", "output");
	cfg.writeFrequency = YAMLReader::optional<Index>(node, "write_frequency", 1);
	cfg.prefix = YAMLReader::optional<std::string>(node, "prefix", "solution");

	const std::string formatStr = YAMLReader::optional<std::string>(node, "format", "vtk");
	if (formatStr == "vtk") {
		cfg.format = residuum::solver::config::OutputConfig::Format::VTK;
	} else if (formatStr == "vtu") {
		cfg.format = residuum::solver::config::OutputConfig::Format::VTU;
	} else {
		throw std::runtime_error("OutputConfigParser: unknown 'format' value '" + formatStr + "'. Valid options: vtk, vtu");
	}

	return cfg;

}
