#ifndef RESIDUUM_SOLVER_STAGE_TRANSIENTCAPABLESTAGE_HPP
#define RESIDUUM_SOLVER_STAGE_TRANSIENTCAPABLESTAGE_HPP

#include <concepts>

#include "core/Types.hpp"
#include "solver/stage/Stage.hpp"

namespace residuum {
	namespace solver {
		namespace stage {

			template<typename ST>
			concept TransientCapableStage = Stage<ST> && requires(ST& stage, Real t, Index step) {

				{ stage.setDt(t) };
				{ stage.setTime(t) };
				{ stage.advance() };
				{ stage.onStepComplete(step, t) };

			}; // concept TransientCapableStage

		} // namespace stage
	} // namespace solver
} // namespace residuum

#endif
