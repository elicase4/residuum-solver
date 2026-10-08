#ifndef RESIDUUM_FEM_EVAL_FIELDSTATE_HPP
#define RESIDUUM_FEM_EVAL_FIELDSTATE_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace residuum {
	namespace fem {
		namespace evaluator {

			template<Index NumDOFs, Index Order>
			struct FieldState {

				Real valueData[NumDOFs];
				Real derivatives[Order == 0 ? 1 : Order][NumDOFs];

				PDE_HOST PDE_DEVICE Real& value(Index dof = 0) { return valueData[dof]; }
				PDE_HOST PDE_DEVICE const Real& value(Index dof = 0) const { return valueData[dof]; }

				PDE_HOST PDE_DEVICE Real& rate(Index dof = 0) requires (Order >= 1) { return derivatives[0][dof]; }
				PDE_HOST PDE_DEVICE const Real& rate(Index dof = 0) const requires (Order >= 1) { return derivatives[0][dof]; }

				PDE_HOST PDE_DEVICE Real& acceleration(Index dof = 0) requires (Order >= 2) { return derivatives[1][dof]; }
				PDE_HOST PDE_DEVICE const Real& acceleration(Index dof = 0) const requires (Order >= 2) { return derivatives[1][dof]; }

			}; // struct FieldState

		} // namespace evaluator
	} // namespace fem
} // namespace residuum

#endif
