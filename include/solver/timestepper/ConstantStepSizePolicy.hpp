#ifndef RESIDUUM_SOLVER_TIMESTEPPER_CONSTANTSTEPSIZEPOLICY_HPP
#define RESIDUUM_SOLVER_TIMESTEPPER_CONSTANTSTEPSIZEPOLICY_HPP

#include "core/Types.hpp"

#include "solver/timestepper/TimeStepSizePolicy.hpp"

namespace residuum {
	namespace solver {
		namespace timestepper {

			class ConstantStepSizePolicy : public TimeStepSizePolicy {
			public:

				explicit ConstantStepSizePolicy(Real dt) : dt_(dt) {}

				Real dt() const override { return dt_; }

				void onStepComplete(Real /*residualNorm*/) override {}

				bool rejectLastStep() const override { return false; }

			private:

				Real dt_;

			}; // class ConstantStepSizePolicy

		} // namespace timestepper
	} // namespace solver
} // namespace residuum

#endif
