#ifndef RESIDUUM_SOLVER_TIMESTEPPER_TIMESTEPSIZEPOLICY_HPP
#define RESIDUUM_SOLVER_TIMESTEPPER_TIMESTEPSIZEPOLICY_HPP

#include "core/Types.hpp"

namespace residuum {
	namespace solver {
		namespace timestepper {

			class TimeStepSizePolicy {
			public:

				virtual ~TimeStepSizePolicy() = default;

				// Size of the step about to be taken.
				virtual Real dt() const = 0;

				// Called once the step just attempted has finished
				virtual void onStepComplete(Real residualNorm) = 0;

				virtual bool rejectLastStep() const = 0;

			}; // class TimeStepSizePolicy

		} // namespace timestepper
	} // namespace solver
} // namespace residuum

#endif
