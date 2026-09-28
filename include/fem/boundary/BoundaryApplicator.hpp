#ifndef RESIDUUM_FEM_BOUNDARY_BOUNDARYAPPLICATOR_HPP
#define RESIDUUM_FEM_BOUNDARY_BOUNDARYAPPLICATOR_HPP

#include <cstring>
#include <unordered_set>

#include "fem/assembly/ElementMap.hpp"
#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/boundary/NaturalBoundaryRegistry.hpp"

#include "fem/evaluator/EvalQuadraturePointBoundary.hpp"
#include "fem/evaluator/EvalQuadraturePointVolume.hpp"

#include "mesh/Mesh.hpp"

#include "topology/TopologicalDOF.hpp"

#include "linalg/types/Vector.hpp"
#include "linalg/types/CSRMatrix.hpp"

namespace residuum {
	namespace fem {
		namespace boundary {

			namespace evaluator = residuum::fem::evaluator;

			template<typename BackendT>
			class BoundaryApplicator {
			public:
				
				template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FormsT, typename QuadratureT>
				requires evaluator::EvalModel<ModelT, EvalQPT>
				static void applyEssentialBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const FormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, linalg::types::Vector<Real, linalg::types::backend::CPU>& F);

				template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointBoundary EvalQPT, typename QuadratureT>
				static void applyNaturalBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const NaturalBoundaryRegistry<EvalQPT>& bcRegistry, const Real time, const EvalEleT& evalEle, const QuadratureT& quadrature, linalg::types::Vector<Real, BackendT>& F);

			}; // class BoundaryApplicator

		} // namespace boundary
	} // namespace fem
} // namespace residuum

#include "backend/cpu/BoundaryApplicator.tpp"
//#include "backend/cuda/BoundaryApplicator.tpp"

#endif
