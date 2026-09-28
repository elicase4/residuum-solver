#include "application/heateq/parser/InitialConditionConfigParser.hpp"

#include <stdexcept>

#include "solver/parser/NodalFieldReadConfigParser.hpp"

residuum::application::heateq::config::InitialConditionConfig residuum::application::heateq::parser::InitialConditionConfigParser::parse(const YAML::Node& node) {

	residuum::application::heateq::config::InitialConditionConfig cfg;

	const YAML::Node& readNode = node["read"];
	if (!readNode) {
		throw std::runtime_error("InitialConditionConfigReader: missing required 'read' section in 'initial_condition'");
	}

	cfg.read = residuum::solver::parser::NodalFieldReadConfigParser::parse(readNode);

	return cfg;

}
