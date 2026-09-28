#ifndef RESIDUUM_FEM_ASSEMBLY_ELEMENTMAP_HPP
#define RESIDUUM_FEM_ASSEMBLY_ELEMENTMAP_HPP

#include "config/Platform.hpp"
#include "core/Types.hpp"

#include "fem/boundary/EssentialBoundaryRegistry.hpp"

#include "topology/TopologicalDOF.hpp"

namespace residuum {
	namespace fem {
		namespace assembly {

			enum class GatherMode { Free, Full, Constrained };
			enum class ScatterMode { Free };

			struct NoSolution {};

			template<Index numDOFs, Index SpatialDim, GatherMode Mode, typename VectorT = NoSolution>
			PDE_HOST PDE_DEVICE void gatherElementVector(const Index* nodeIDs, Index nodesPerElement, const Real* nodeCoords, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry* bcRegistry, const Real time, Real* Ue, const VectorT* U = nullptr);

			template<Index numDOFs, ScatterMode Mode, typename VectorT>
			PDE_HOST PDE_DEVICE void scatterElementVector(const Index* nodeIDs, Index nodesPerElement, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real* Xe, VectorT& X, Real coefficient = Real(1));

			// scatters only Ke's diagonal entries (Ke[a,a]) into a plain vector -- no sparsity
			// structure involved, so unlike scatterElementMatrix this needs no merge/search at all
			template<Index numDOFs, ScatterMode Mode, typename VectorT>
			PDE_HOST PDE_DEVICE void scatterElementDiagonal(const Index* nodeIDs, Index nodesPerElement, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real* Ke, VectorT& diag);

			template<Index numDOFs, Index MaxNodesPerElement, ScatterMode Mode, typename MatrixT>
			PDE_HOST PDE_DEVICE void scatterElementMatrix(const Index* nodeIDs, Index nodesPerElement, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real* Ke, MatrixT& K, Real coefficient = Real(1));

		} // namespace assembly
	} // namespace fem
} // namespace residuum

#include "fem/assembly/backend/cpu/ElementMap.tpp"

#endif
