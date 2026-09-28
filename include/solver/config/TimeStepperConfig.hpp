#ifndef RESIDUUM_SOLVER_CONFIG_TIMESTEPPERCONFIG_HPP
#define RESIDUUM_SOLVER_CONFIG_TIMESTEPPERCONFIG_HPP

#include "core/Types.hpp"
#include "solver/config/NonlinearSolverConfig.hpp"
#include "solver/config/LinearSolverConfig.hpp"
#include "solver/config/TimeStepSizeConfig.hpp"

namespace residuum {
	namespace solver {
		namespace config {
			
			struct TimeStepperConfig {

				enum class Type {
					ForwardEuler,
					BackwardEuler,
					GeneralizedAlpha,
					RK4
				}; // enum class Type
	
				Type type = Type::ForwardEuler;
				
				// Time interval
				Real t0 = 0.0;
				Real tf = 1.0;

				TimeStepSizeConfig stepSize;

				Real rhoInf = 0.5;

				// Solver config
				NonlinearSolverConfig nonlinearSolver;
				LinearSolverConfig linearSolver;

			}; // struct TimeStepperConfig

		} // namespace config
	} // namespace solver
} // namespace residuum

#endif
