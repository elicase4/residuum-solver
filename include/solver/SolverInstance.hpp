#ifndef RESIDUUM_SOLVER_SOLVERINSTANCE_HPP
#define RESIDUUM_SOLVER_SOLVERINSTANCE_HPP

#include "solver/config/SolverConfig.hpp"

namespace residuum {
	namespace solver {

		enum class SolverMode {
			SteadyLinear,
			SteadyNonlinear,
			TransientLinear,
			TransientNonlinear
		};

		struct SolverInstance {

			SolverMode mode;

			const config::LinearSolverConfig* linear = nullptr;

			const config::NonlinearSolverConfig* nonlinear = nullptr;

			const config::TimeStepperConfig* timestepper = nullptr;

		}; // struct SolverInstance

		SolverInstance resolveSolverInstance(const config::SolverConfig& solver);

		inline bool isNonlinear(SolverMode mode) {
			return mode == SolverMode::SteadyNonlinear || mode == SolverMode::TransientNonlinear;
		}

		inline bool isTransient(SolverMode mode) {
			return mode == SolverMode::TransientLinear || mode == SolverMode::TransientNonlinear;
		}

	} // namespace solver
} // namespace residuum

#endif
