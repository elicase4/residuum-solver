#ifndef RESIDUUM_FEM_EVAL_STATEVECTOR_HPP
#define RESIDUUM_FEM_EVAL_STATEVECTOR_HPP

#include <concepts>

#include "core/Types.hpp"
#include "topology/TopologicalDOF.hpp"

namespace residuum {
	namespace fem {
		namespace evaluator {

			// constrainedValues is the flat, per-call-lowered Dirichlet buffer (see fem::boundary::lowerEssentialBoundaryValues), sized topoDOF.numGlobalDOFs()
			template<typename StateVectorT, Index numDOFs, typename QuadraturePointT>
			concept StateVector = requires (StateVectorT& sv, QuadraturePointT& qp,
				const Index* nodeIDs, Index nodesPerElement, const Real* nodeCoords,
				const topology::TopologicalDOF<numDOFs>& topoDOF,
				const Real* constrainedValues, Real time) {

				{ sv.gatherElement(nodeIDs, nodesPerElement, nodeCoords, topoDOF, constrainedValues, time) } -> std::same_as<void>;
				{ sv.interpolate(qp) } -> std::same_as<void>;

			}; // concept StateVector

		} // namespace evaluator
	} // namespace fem
} // namespace residuum

#endif
