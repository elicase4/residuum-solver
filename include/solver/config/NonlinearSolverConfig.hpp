#ifndef RESIDUUM_SOLVER_CONFIG_NONLINEARSOLVERCONFIG_HPP
#define RESIDUUM_SOLVER_CONFIG_NONLINEARSOLVERCONFIG_HPP

#include "core/Types.hpp"
#include "solver/config/LinearSolverConfig.hpp"

namespace residuum {
	namespace solver {
		namespace config {

			struct NonlinearSolverConfig {

				enum class Type {
					Newton,
					Picard
				};

				Type type = Type::Newton;

				Real absoluteTolerance = 1e-10;

				Real relativeTolerance = 1e-8;

				Index maxIterations = 50;

				LinearSolverConfig linearSolver;

			}; // struct NonlinearSolverConfig

		} // namespace config
	} // namespace solver
} // namespace residuum

#endif
