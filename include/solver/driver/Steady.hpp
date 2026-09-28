#ifndef RESIDUUM_SOLVER_DRIVER_STEADY_HPP
#define RESIDUUM_SOLVER_DRIVER_STEADY_HPP

#include "solver/driver/Driver.hpp"
#include "solver/stage/Stage.hpp"

namespace residuum {
	namespace solver {
		namespace driver {

			template<stage::Stage StageT>
			class Steady {
			public:

				Steady() {
					static_assert(Driver<Steady, StageT>);
				}

				bool solve(StageT& stage) {

					stage.initialize();
					stage.assemble();
					bool converged = stage.solve();
					stage.finalize();

					return converged;
				}

			}; // class Steady

		} // namespace driver
	} // namespace solver
} // namespace residuum

#endif
