#include "application/heateq/parser/OutputFieldConfigParser.hpp"

#include <stdexcept>

residuum::application::heateq::config::OutputFieldQuantity residuum::application::heateq::parser::OutputFieldConfigParser::parseQuantity(const std::string& str) {

	if (str == "heat_flux") return residuum::application::heateq::config::OutputFieldQuantity::HeatFlux;
	if (str == "conductivity") return residuum::application::heateq::config::OutputFieldQuantity::Conductivity;
	if (str == "specific_heat") return residuum::application::heateq::config::OutputFieldQuantity::SpecificHeat;
	if (str == "residual") return residuum::application::heateq::config::OutputFieldQuantity::Residual;
	throw std::runtime_error("OutputFieldConfigParser: unknown output field '" + str + "'");

}

residuum::application::heateq::config::OutputFieldConfig residuum::application::heateq::parser::OutputFieldConfigParser::parse(const YAML::Node& node) {

	residuum::application::heateq::config::OutputFieldConfig cfg;

	for (const auto& fieldNode : node) {
		cfg.fields.push_back(residuum::application::heateq::parser::OutputFieldConfigParser::parseQuantity(fieldNode.as<std::string>()));
	}

	return cfg;

}
