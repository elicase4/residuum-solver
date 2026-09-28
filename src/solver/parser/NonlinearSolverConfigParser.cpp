#include "solver/parser/NonlinearSolverConfigParser.hpp"
#include "solver/parser/LinearSolverConfigParser.hpp"
#include "io/YAMLReader.hpp"

residuum::solver::config::NonlinearSolverConfig::Type residuum::solver::parser::NonlinearSolverConfigParser::parseNonlinearSolverType(const std::string& str) {

	if (str == "newton")
		return residuum::solver::config::NonlinearSolverConfig::Type::Newton;
	if (str == "picard")
		return residuum::solver::config::NonlinearSolverConfig::Type::Picard;

	throw std::runtime_error("Unknown nonlinear solver type: '" + str + "'. Valid options: newton, picard");

}

residuum::solver::config::NonlinearSolverConfig residuum::solver::parser::NonlinearSolverConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	residuum::solver::config::NonlinearSolverConfig cfg;

	cfg.type = parseNonlinearSolverType(YAMLReader::required<std::string>(node, "type"));
	cfg.absoluteTolerance = YAMLReader::optional<Real>(node, "absolute_tolerance", 1e-10);
	cfg.relativeTolerance = YAMLReader::optional<Real>(node, "relative_tolerance", 1e-8);
	cfg.maxIterations = YAMLReader::optional<Index>(node, "max_iterations", 50);

	if (node["linear_solver"]) {
		cfg.linearSolver = LinearSolverConfigParser::parse(node["linear_solver"]);
	}

	return cfg;
}
