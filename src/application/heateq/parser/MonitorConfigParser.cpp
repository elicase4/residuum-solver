#include "application/heateq/parser/MonitorConfigParser.hpp"

#include <stdexcept>

#include "io/YAMLReader.hpp"

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

	const std::string reductionStr = YAMLReader::required<std::string>(node, "reduction");
	if (reductionStr == "integral") {
		cfg.reduction = residuum::fem::quantity::Reduction::Integral;
	} else if (reductionStr == "average") {
		cfg.reduction = residuum::fem::quantity::Reduction::Average;
	} else {
		throw std::runtime_error("MonitorConfigParser: monitor '" + cfg.name + "': unknown reduction '" + reductionStr + "'");
	}

	const YAML::Node& termsNode = node["terms"];
	if (!termsNode) {
		throw std::runtime_error("MonitorConfigParser: monitor '" + cfg.name + "': missing required 'terms' list");
	}

	for (const auto& termNode : termsNode) {

		residuum::application::heateq::config::MonitorTermConfig term;
		term.boundary = YAMLReader::required<Int>(termNode, "boundary");
		term.coefficient = YAMLReader::optional<Real>(termNode, "coefficient", Real(1.0));
		cfg.terms.push_back(term);

	}

	const YAML::Node& outputNode = node["output"];
	if (outputNode) {
		cfg.output.console = YAMLReader::optional<bool>(outputNode, "console", true);
		cfg.output.file = YAMLReader::optional<std::string>(outputNode, "file", "");
	}

	return cfg;

}
