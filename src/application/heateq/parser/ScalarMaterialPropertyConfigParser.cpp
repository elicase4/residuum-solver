#include "application/heateq/parser/ScalarMaterialPropertyConfigParser.hpp"

#include "io/YAMLReader.hpp"

residuum::application::heateq::config::ScalarMaterialPropertyConfig residuum::application::heateq::parser::ScalarMaterialPropertyConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	residuum::application::heateq::config::ScalarMaterialPropertyConfig cfg;

	cfg.value = YAMLReader::required<Real>(node, "value");
	cfg.unit = YAMLReader::required<std::string>(node, "unit");

	return cfg;

}
