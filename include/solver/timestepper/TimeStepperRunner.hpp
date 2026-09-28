#ifndef RESIDUUM_SOLVER_TIMESTEPPER_TIMESTEPPERRUNNER_HPP
#define RESIDUUM_SOLVER_TIMESTEPPER_TIMESTEPPERRUNNER_HPP

#include "core/Types.hpp"

namespace residuum {
	namespace solver {
		namespace timestepper {

			class TimeStepperRunner {
			public:

				virtual ~TimeStepperRunner() = default;

				// Advance one timestep. Returns false if that step's solve failed.
				virtual bool step() = 0;

				virtual bool finished() const = 0;

				virtual Index currentStep() const = 0;

				virtual Real currentTime() const = 0;

			}; // class TimeStepperRunner

		} // namespace timestepper
	} // namespace solver
} // namespace residuum

#endif
