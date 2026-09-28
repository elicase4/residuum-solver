#include "solver/parser/MeshConfigParser.hpp"
#include "io/YAMLReader.hpp"

residuum::solver::config::MeshConfig residuum::solver::parser::MeshConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	residuum::solver::config::MeshConfig cfg;

	cfg.file = YAMLReader::required<std::string>(node, "file");

	return cfg;

}
