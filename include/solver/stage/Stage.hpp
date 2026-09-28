#ifndef RESIDUUM_SOLVER_STAGE_STAGE_HPP
#define RESIDUUM_SOLVER_STAGE_STAGE_HPP

#include <concepts>

namespace residuum {
	namespace solver {
		namespace stage {

			template<typename ST>
			concept Stage = requires(ST& stage) {

				{ stage.initialize() };
				{ stage.assemble() };
				{ stage.solve() } -> std::same_as<bool>;
				{ stage.finalize() };

			}; // concept Stage

		} // namespace stage
	} // namespace solver
} // namespace residuum

#endif
