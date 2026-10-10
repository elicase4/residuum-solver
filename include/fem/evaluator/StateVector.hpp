#ifndef RESIDUUM_FEM_EVAL_STATEVECTOR_HPP
#define RESIDUUM_FEM_EVAL_STATEVECTOR_HPP

#include <concepts>

#include "core/Types.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "topology/TopologicalDOF.hpp"

namespace residuum {
	namespace fem {
		namespace evaluator {

			template<typename StateVectorT, Index numDOFs, typename QuadraturePointT>
			concept StateVector = requires (StateVectorT& sv, QuadraturePointT& qp,
				const Index* nodeIDs, Index nodesPerElement, const Real* nodeCoords,
				const topology::TopologicalDOF<numDOFs>& topoDOF,
				const fem::boundary::EssentialBoundaryRegistry* bcRegistry, Real time) {

				{ sv.gatherElement(nodeIDs, nodesPerElement, nodeCoords, topoDOF, bcRegistry, time) } -> std::same_as<void>;
				{ sv.interpolate(qp) } -> std::same_as<void>;

			}; // concept StateVector

		} // namespace evaluator
	} // namespace fem
} // namespace residuum

#endif
