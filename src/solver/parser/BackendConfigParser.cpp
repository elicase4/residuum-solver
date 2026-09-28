#include "solver/parser/BackendConfigParser.hpp"
#include "io/YAMLReader.hpp"

#include <stdexcept>

residuum::solver::config::BackendConfig::Type residuum::solver::parser::BackendConfigParser::parseBackendType(const std::string& str) {

	if (str == "cpu")
		return config::BackendConfig::Type::CPU;
	if (str == "cuda")
		return config::BackendConfig::Type::CUDA;

	throw std::runtime_error("Unknown backend type: '" + str + "'. Valid options: cpu, cuda");

}

residuum::solver::config::BackendConfig residuum::solver::parser::BackendConfigParser::parse(const YAML::Node& root) {

	using io::YAMLReader;

	config::BackendConfig cfg;

	cfg.type = parseBackendType(YAMLReader::optional<std::string>(root, "backend", std::string("cpu")));

	return cfg;

}
