#include "solver/SolverInstance.hpp"

#include <stdexcept>

residuum::solver::SolverInstance residuum::solver::resolveSolverInstance(const residuum::solver::config::SolverConfig& solver) {

	SolverInstance instance;

	const bool transient = (solver.driver.type == config::DriverConfig::Type::Transient);

	if (transient && !solver.timestepper.has_value()) {
		throw std::runtime_error("SolverInstance: transient driver requires a 'timestepper' section");
	}

	if (solver.nonlinear.has_value()) {

		instance.mode = transient ? SolverMode::TransientNonlinear : SolverMode::SteadyNonlinear;
		instance.nonlinear = &(*solver.nonlinear);
		instance.linear = &solver.nonlinear->linearSolver;

	} else if (solver.linear.has_value()) {

		instance.mode = transient ? SolverMode::TransientLinear : SolverMode::SteadyLinear;
		instance.linear = &(*solver.linear);

	} else {
		throw std::runtime_error("SolverInstance: solver config requires either 'linear' or 'nonlinear'");
	}

	if (transient) {
		instance.timestepper = &(*solver.timestepper);
	}

	return instance;

}
