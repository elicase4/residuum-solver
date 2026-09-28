#ifndef RESIDUUM_SOLVER_DRIVER_DRIVER_HPP
#define RESIDUUM_SOLVER_DRIVER_DRIVER_HPP

#include <concepts>

namespace residuum {
	namespace solver {
		namespace driver {

			template<typename DT, typename StageT>
			concept Driver = requires(DT& driver, StageT& stage) {
				
				{ driver.solve(stage) } -> std::same_as<bool>;

			}; // concept Driver

		} // namespace driver
	} // namespace solver
} // namespace residuum

#endif
