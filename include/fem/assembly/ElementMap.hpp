#ifndef RESIDUUM_FEM_ASSEMBLY_ELEMENTMAP_HPP
#define RESIDUUM_FEM_ASSEMBLY_ELEMENTMAP_HPP

#include "config/Platform.hpp"
#include "core/Types.hpp"

#include "mesh/Mesh.hpp"

#include "topology/TopologicalDOF.hpp"

namespace residuum {
	namespace fem {
		namespace assembly {

			enum class GatherMode { Free, Full, Constrained };
			enum class ScatterMode { Free };

			struct EmptySolution {};

			// extracts this element's nodes' physical coordinates into a flat, SpatialDim-major buffer, factoring out the loop every assembly/evaluation element pass repeated
			template<Index SpatialDim>
			PDE_HOST PDE_DEVICE void gatherElementCoords(const mesh::Mesh& mesh, const Index* nodeIDs, Index nodesPerElement, Real* nodeCoords);

			// materializes a quadrature rule's points/weights into flat buffers, factoring out the pair every assembly/evaluation element pass repeated
			template<typename QuadratureT>
			PDE_HOST PDE_DEVICE void materializeQuadrature(const QuadratureT& quadrature, Real* xi, Real* w) {
				quadrature.getPoints(xi);
				quadrature.getWeights(w);
			}

			// gathers one element's Ue buffer; constrainedValues (sized topoDOF.numGlobalDOFs(), see fem::boundary::lowerEssentialBoundaryValues) supplies constrained-DOF values for Full/Constrained modes -- Free mode never reads it
			template<Index numDOFs, GatherMode Mode, typename VectorT = EmptySolution>
			PDE_HOST PDE_DEVICE void gatherElementVector(const Index* nodeIDs, Index nodesPerElement, const topology::TopologicalDOF<numDOFs>& topoDOF, Real* Ue, const VectorT* U = nullptr, const Real* constrainedValues = nullptr);

			template<Index numDOFs, ScatterMode Mode, typename VectorT>
			PDE_HOST PDE_DEVICE void scatterElementVector(const Index* nodeIDs, Index nodesPerElement, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real* Xe, VectorT& X, Real coefficient = Real(1));

			// scatters only Ke's diagonal entries (Ke[a,a]) into a plain vector, so unlike scatterElementMatrix this needs no sparsity merge/search at all
			template<Index numDOFs, ScatterMode Mode, typename VectorT>
			PDE_HOST PDE_DEVICE void scatterElementDiagonal(const Index* nodeIDs, Index nodesPerElement, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real* Ke, VectorT& diag);

			template<Index numDOFs, Index MaxNodesPerElement, ScatterMode Mode, typename MatrixT>
			PDE_HOST PDE_DEVICE void scatterElementMatrix(const Index* nodeIDs, Index nodesPerElement, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real* Ke, MatrixT& K, Real coefficient = Real(1));

		} // namespace assembly
	} // namespace fem
} // namespace residuum

#include "fem/assembly/backend/cpu/ElementMap.tpp"

#endif
