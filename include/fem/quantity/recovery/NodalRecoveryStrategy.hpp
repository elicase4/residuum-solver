#ifndef RESIDUUM_FEM_QUANTITY_RECOVERY_NODALRECOVERYSTRATEGY_HPP
#define RESIDUUM_FEM_QUANTITY_RECOVERY_NODALRECOVERYSTRATEGY_HPP

#include "core/Types.hpp"

namespace residuum {
	namespace fem {
		namespace quantity {
			namespace recovery {

				// contribution/nodalAccum both have length numComponents; nodalWeight is a single Real
				template<typename StrategyT>
				concept NodalRecoveryStrategy = requires (StrategyT s, const Real* contribution, Real weight, Real* nodalAccum, Real* nodalWeight, Real nodalWeightValue, Index numComponents) {
					{ s.accumulate(contribution, weight, nodalAccum, nodalWeight, numComponents) } -> std::same_as<void>;
					{ s.finalize(nodalAccum, nodalWeightValue, numComponents) } -> std::same_as<void>;
				}; // concept NodalRecoveryStrategy

			} // namespace recovery
		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
