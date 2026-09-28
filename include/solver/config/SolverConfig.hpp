#ifndef RESIDUUM_SOLVER_CONFIG_SOLVERCONFIG_HPP
#define RESIDUUM_SOLVER_CONFIG_SOLVERCONFIG_HPP

#include <optional>

#include "solver/config/DriverConfig.hpp"
#include "solver/config/TimeStepperConfig.hpp"
#include "solver/config/NonlinearSolverConfig.hpp"
#include "solver/config/LinearSolverConfig.hpp"

namespace residuum {
	namespace solver {
		namespace config {

			struct SolverConfig {
				
				DriverConfig driver;
				
				std::optional<TimeStepperConfig> timestepper;

				std::optional<NonlinearSolverConfig> nonlinear;

				std::optional<LinearSolverConfig> linear;

			}; // struct SolverConfig

		} // namespace config
	} // namespace solver
} // namespace residuum

#endif
