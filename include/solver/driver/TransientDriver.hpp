#ifndef RESIDUUM_SOLVER_DRIVER_TRANSIENTDRIVER_HPP
#define RESIDUUM_SOLVER_DRIVER_TRANSIENTDRIVER_HPP

#include <concepts>

#include "solver/timestepper/TimeStepperRunner.hpp"

namespace residuum {
	namespace solver {
		namespace driver {

			template<typename DT, typename StageT>
			concept TransientDriver = requires(DT& driver, StageT& stage, timestepper::TimeStepperRunner& stepper) {

				{ driver.solve(stage, stepper) } -> std::same_as<bool>;

			}; // concept TransientDriver

		} // namespace driver
	} // namespace solver
} // namespace residuum

#endif
