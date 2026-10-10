#include "application/heateq/parser/MonitorConfigParser.hpp"

#include <stdexcept>

#include "io/YAMLReader.hpp"
#include "solver/parser/EvaluationMethodParser.hpp"
#include "solver/parser/MonitorOutputConfigParser.hpp"
#include "solver/parser/MonitorTermConfigParser.hpp"
#include "solver/parser/ReductionParser.hpp"

residuum::application::heateq::config::MonitorConfig::Quantity residuum::application::heateq::parser::MonitorConfigParser::parseQuantity(const std::string& str) {

	if (str == "heat_flux"){
		return residuum::application::heateq::config::MonitorConfig::Quantity::HeatFlux;
	} else {
		throw std::runtime_error("Unknown monitor quantity: " + str);
	}

}

residuum::application::heateq::config::MonitorConfig residuum::application::heateq::parser::MonitorConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	residuum::application::heateq::config::MonitorConfig cfg;

	cfg.name = YAMLReader::required<std::string>(node, "name");
	cfg.quantity = residuum::application::heateq::parser::MonitorConfigParser::parseQuantity(YAMLReader::required<std::string>(node, "quantity"));
	cfg.fluxMethod = residuum::solver::parser::EvaluationMethodParser::parse(YAMLReader::optional<std::string>(node, "flux_method", "weak"));
	cfg.reduction = residuum::solver::parser::ReductionParser::parse(YAMLReader::required<std::string>(node, "reduction"));

	const YAML::Node& termsNode = node["terms"];
	if (!termsNode) {
		throw std::runtime_error("MonitorConfigParser: monitor '" + cfg.name + "': missing required 'terms' list");
	}

	for (const auto& termNode : termsNode) {
		cfg.terms.push_back(residuum::solver::parser::MonitorTermConfigParser::parse(termNode));
	}

	const YAML::Node& outputNode = node["output"];
	if (outputNode) {
		cfg.output = residuum::solver::parser::MonitorOutputConfigParser::parse(outputNode);
	}

	return cfg;

}
