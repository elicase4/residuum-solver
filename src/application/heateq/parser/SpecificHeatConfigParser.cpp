#include "application/heateq/parser/SpecificHeatConfigParser.hpp"
#include "io/YAMLReader.hpp"

residuum::application::heateq::config::SpecificHeatConfig::Type residuum::application::heateq::parser::SpecificHeatConfigParser::parseSpecificHeatType(const std::string& str) {

	if (str == "constant") {
		return residuum::application::heateq::config::SpecificHeatConfig::Type::Constant;
	}

	if (str == "temperature_dependent") {
		return residuum::application::heateq::config::SpecificHeatConfig::Type::TemperatureDependent;
	}

	throw std::runtime_error("Unknown specific heat type: " + str);

}

residuum::application::heateq::config::SpecificHeatConfig residuum::application::heateq::parser::SpecificHeatConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;
	using Type = residuum::application::heateq::config::SpecificHeatConfig::Type;

	residuum::application::heateq::config::SpecificHeatConfig cfg;

	cfg.type = residuum::application::heateq::parser::SpecificHeatConfigParser::parseSpecificHeatType(YAMLReader::optional<std::string>(node, "type", "constant"));

	if (cfg.type == Type::Constant) {
		cfg.value = YAMLReader::required<Real>(node, "value");
	} else if (cfg.type == Type::TemperatureDependent) {
		cfg.valueExpression = YAMLReader::required<std::string>(node, "value_expression");
		cfg.gradientExpression = YAMLReader::required<std::string>(node, "gradient_expression");
	}

	cfg.unit = YAMLReader::required<std::string>(node, "unit");

	return cfg;

}
