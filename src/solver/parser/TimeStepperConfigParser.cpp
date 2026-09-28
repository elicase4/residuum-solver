#include "solver/parser/TimeStepperConfigParser.hpp"
#include "solver/parser/NonlinearSolverConfigParser.hpp"
#include "solver/parser/LinearSolverConfigParser.hpp"
#include "solver/parser/TimeStepSizeConfigParser.hpp"
#include "io/YAMLReader.hpp"

residuum::solver::config::TimeStepperConfig::Type residuum::solver::parser::TimeStepperConfigParser::parseTimeStepperType(const std::string& str) {

	if (str == "forward_euler")
		return residuum::solver::config::TimeStepperConfig::Type::ForwardEuler;
	if (str == "backward_euler")
		return residuum::solver::config::TimeStepperConfig::Type::BackwardEuler;
	if (str == "generalized_alpha") 
		return residuum::solver::config::TimeStepperConfig::Type::GeneralizedAlpha;
	if (str == "rk4")
		return residuum::solver::config::TimeStepperConfig::Type::RK4;

	throw std::runtime_error("Unknown time stepper type: '" + str + "'. Valid options: forward_euler, backward_euler, generalized_alpha, rk4");

}

residuum::solver::config::TimeStepperConfig residuum::solver::parser::TimeStepperConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	residuum::solver::config::TimeStepperConfig cfg;

	cfg.type = parseTimeStepperType(YAMLReader::required<std::string>(node, "type"));
	cfg.t0 = YAMLReader::optional<Real>(node, "t0", 0.0);
	cfg.tf = YAMLReader::optional<Real>(node, "tf", 0.0);
	cfg.rhoInf = YAMLReader::optional<Real>(node, "rho_inf", 0.5);

	const YAML::Node& stepSize = node["step_size"];
	if (!stepSize) {
		throw std::runtime_error("TimeStepperConfigParser: missing required 'step_size' section");
	}
	cfg.stepSize = TimeStepSizeConfigParser::parse(stepSize);

	if (node["nonlinear_solver"]) {
		cfg.nonlinearSolver = NonlinearSolverConfigParser::parse(node["nonlinear_solver"]);
	}

	if (node["linear_solver"]) {
		cfg.linearSolver = LinearSolverConfigParser::parse(node["linear_solver"]);
	}

	return cfg;

}
