#include "solver/parser/SolverConfigParser.hpp"
#include "solver/parser/DriverConfigParser.hpp"
#include "solver/parser/TimeStepperConfigParser.hpp"
#include "solver/parser/LinearSolverConfigParser.hpp"
#include "solver/parser/NonlinearSolverConfigParser.hpp"
#include "io/YAMLReader.hpp"

residuum::solver::config::SolverConfig residuum::solver::parser::SolverConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	residuum::solver::config::SolverConfig cfg;

	const YAML::Node& driver = node["driver"];
	if (!driver) {
		throw std::runtime_error("SolverConfigParser: missing required 'driver' section");
	}
	cfg.driver = residuum::solver::parser::DriverConfigParser::parse(driver);

	const YAML::Node& timestepper = node["timestepper"];
	if (timestepper) {
		cfg.timestepper = residuum::solver::parser::TimeStepperConfigParser::parse(timestepper);
	}

	const YAML::Node& nonlinear = node["nonlinear"];
	if (nonlinear) {
		cfg.nonlinear = residuum::solver::parser::NonlinearSolverConfigParser::parse(nonlinear);
	}

	const YAML::Node& linear = node["linear"];
	if (linear) {
		cfg.linear = residuum::solver::parser::LinearSolverConfigParser::parse(linear);
	}

	return cfg;

}
