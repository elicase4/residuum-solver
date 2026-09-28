#ifndef RESIDUUM_SOLVER_DRIVER_TRANSIENT_HPP
#define RESIDUUM_SOLVER_DRIVER_TRANSIENT_HPP

#include "solver/driver/TransientDriver.hpp"
#include "solver/stage/Stage.hpp"
#include "solver/timestepper/TimeStepperRunner.hpp"

namespace residuum {
	namespace solver {
		namespace driver {

			// owns only the outer loop; each step's assemble/solve/advance/output lives in the timestepper
			template<stage::Stage StageT>
			class Transient {
			public:

				Transient() {
					static_assert(TransientDriver<Transient, StageT>);
				}

				bool solve(StageT& stage, timestepper::TimeStepperRunner& stepper) {

					stage.initialize();

					while (!stepper.finished()) {
						if (!stepper.step()) {
							stage.finalize();
							return false;
						}
					}

					stage.finalize();
					return true;

				}

			}; // class Transient

		} // namespace driver
	} // namespace solver
} // namespace residuum

#endif
