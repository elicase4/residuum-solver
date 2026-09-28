#include "application/heateq/parser/BoundaryConditionConfigParser.hpp"
#include "application/heateq/parser/ConductivityConfigParser.hpp"
#include "io/YAMLReader.hpp"
#include "solver/parser/NodalFieldReadConfigParser.hpp"

residuum::application::heateq::config::BoundaryConditionConfig::Type residuum::application::heateq::parser::BoundaryConditionConfigParser::parseBoundaryConditionType(const std::string& str) {

	if (str == "value") {
		return residuum::application::heateq::config::BoundaryConditionConfig::Type::Value;
	}

	if (str == "flux") {
		return residuum::application::heateq::config::BoundaryConditionConfig::Type::Flux;
	}

	throw std::runtime_error("Unknown boundary type: " + str);

}

residuum::application::heateq::config::BoundaryConditionConfig::Form residuum::application::heateq::parser::BoundaryConditionConfigParser::parseBoundaryConditionForm(const std::string& str) {

	if (str == "flux_bc") {
		return residuum::application::heateq::config::BoundaryConditionConfig::Form::FluxBC;
	}

	if (str == "value_bc") {
		return residuum::application::heateq::config::BoundaryConditionConfig::Form::ValueBC;
	}

	throw std::runtime_error("Unknown boundary form: " + str);

}

residuum::application::heateq::config::BoundaryConditionConfig residuum::application::heateq::parser::BoundaryConditionConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;
	using residuum::solver::config::NodalFieldReadConfig;
	using residuum::solver::parser::NodalFieldReadConfigParser;

	residuum::application::heateq::config::BoundaryConditionConfig cfg;

	cfg.boundaryID = YAMLReader::required<Int>(node, "boundary");
	cfg.type = BoundaryConditionConfigParser::parseBoundaryConditionType(YAMLReader::required<std::string>(node, "type"));

	const YAML::Node& readNode = node["read"];
	if (!readNode) {
		throw std::runtime_error("BoundaryConditionConfigReader: missing required 'read' section in a boundary_conditions entry");
	}

	// only the Mode enum is reused here; file/expression stay local fields since Flux needs one expression per spatial component
	cfg.mode = NodalFieldReadConfigParser::parseMode(YAMLReader::optional<std::string>(readNode, "mode", "expression"));

	if (cfg.mode == NodalFieldReadConfig::Mode::File) {
		cfg.file = YAMLReader::required<std::string>(readNode, "file");
	} else if (cfg.type == residuum::application::heateq::config::BoundaryConditionConfig::Type::Flux) {
		cfg.fluxExpression = YAMLReader::required<std::vector<std::string>>(readNode, "expression");
	} else {
		cfg.expression = YAMLReader::required<std::string>(readNode, "expression");
	}

	cfg.unit = YAMLReader::required<std::string>(readNode, "unit");

	for (auto& form : YAMLReader::required<std::vector<std::string>>(node, "forms")) {
		cfg.forms.push_back(BoundaryConditionConfigParser::parseBoundaryConditionForm(form));
	}
	cfg.model = residuum::application::heateq::parser::ConductivityConfigParser::parseConductivityType(YAMLReader::required<std::string>(node, "model"));

	return cfg;

}
