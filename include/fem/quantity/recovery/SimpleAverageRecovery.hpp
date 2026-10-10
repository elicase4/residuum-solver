#ifndef RESIDUUM_FEM_QUANTITY_RECOVERY_SIMPLEAVERAGERECOVERY_HPP
#define RESIDUUM_FEM_QUANTITY_RECOVERY_SIMPLEAVERAGERECOVERY_HPP

#include "core/Types.hpp"

namespace residuum {
	namespace fem {
		namespace quantity {
			namespace recovery {

				// volume-and-shape-function-weighted mean of every touching element's quadrature-point contribution at a node (a lumped projection, not a true L2 projection)
				struct SimpleAverageRecovery {

					void accumulate(const Real* contribution, Real weight, Real* nodalAccum, Real* nodalWeight, Index numComponents) const {
						for (Index c = 0; c < numComponents; ++c) nodalAccum[c] += weight * contribution[c];
						*nodalWeight += weight;
					}

					void finalize(Real* nodalAccum, Real nodalWeight, Index numComponents) const {
						for (Index c = 0; c < numComponents; ++c) nodalAccum[c] /= nodalWeight;
					}

				}; // struct SimpleAverageRecovery

			} // namespace recovery
		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
