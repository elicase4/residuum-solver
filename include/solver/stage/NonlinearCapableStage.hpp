#ifndef RESIDUUM_SOLVER_STAGE_NONLINEARCAPABLESTAGE_HPP
#define RESIDUUM_SOLVER_STAGE_NONLINEARCAPABLESTAGE_HPP

#include <concepts>

#include "core/Types.hpp"
#include "solver/stage/Stage.hpp"

namespace residuum {
	namespace solver {
		namespace stage {

			template<typename ST>
			concept NonlinearCapableStage = Stage<ST> && requires(ST& stage) {

				{ stage.residualNorm() } -> std::convertible_to<Real>;

				// the raw residual vector residualNorm() just assembled, for per-DOF reporting
				{ stage.residual() };

				{ stage.solveLinearStep() } -> std::same_as<bool>;

			}; // concept NonlinearCapableStage

		} // namespace stage
	} // namespace solver
} // namespace residuum

#endif
