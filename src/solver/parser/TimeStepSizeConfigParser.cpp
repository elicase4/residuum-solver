#include "solver/parser/TimeStepSizeConfigParser.hpp"

#include <stdexcept>

#include "io/YAMLReader.hpp"

residuum::solver::config::TimeStepSizeConfig::Mode residuum::solver::parser::TimeStepSizeConfigParser::parseMode(const std::string& str) {

	if (str == "constant")
		return config::TimeStepSizeConfig::Mode::Constant;
	if (str == "adaptive")
		return config::TimeStepSizeConfig::Mode::Adaptive;
	if (str == "pseudo_transient")
		return config::TimeStepSizeConfig::Mode::PseudoTransient;

	throw std::runtime_error("Unknown timestepper step_size mode: '" + str + "'. Valid options: constant, adaptive, pseudo_transient");

}

residuum::solver::config::TimeStepSizeConfig residuum::solver::parser::TimeStepSizeConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	config::TimeStepSizeConfig cfg;

	cfg.mode = parseMode(YAMLReader::required<std::string>(node, "mode"));
	cfg.dt = YAMLReader::required<Real>(node, "dt");

	// mode-specific knobs are optional; nothing reads them yet regardless of mode
	if (node["dt_min"]) cfg.dtMin = YAMLReader::required<Real>(node, "dt_min");
	if (node["dt_max"]) cfg.dtMax = YAMLReader::required<Real>(node, "dt_max");
	if (node["growth_factor"]) cfg.growthFactor = YAMLReader::required<Real>(node, "growth_factor");
	if (node["shrink_factor"]) cfg.shrinkFactor = YAMLReader::required<Real>(node, "shrink_factor");
	if (node["error_tolerance"]) cfg.errorTolerance = YAMLReader::required<Real>(node, "error_tolerance");
	if (node["convergence_tolerance"]) cfg.convergenceTolerance = YAMLReader::required<Real>(node, "convergence_tolerance");

	return cfg;

}
