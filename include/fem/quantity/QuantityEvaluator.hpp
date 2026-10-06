#ifndef RESIDUUM_FEM_QUANTITY_QUANTITYEVALUATOR_HPP
#define RESIDUUM_FEM_QUANTITY_QUANTITYEVALUATOR_HPP

#include <array>

#include "config/Platform.hpp"

#include "core/Types.hpp"

#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/evaluator/EvalElement.hpp"
#include "fem/evaluator/EvalModel.hpp"
#include "fem/evaluator/EvalQuadraturePointBoundary.hpp"
#include "fem/evaluator/EvalQuadraturePointVolume.hpp"
#include "fem/quantity/BoundaryQuantityRegistry.hpp"
#include "fem/quantity/Reduction.hpp"
#include "fem/quantity/QuantityForms.hpp"
#include "fem/quantity/recovery/NodalRecoveryStrategy.hpp"

#include "mesh/Mesh.hpp"

#include "linalg/types/Vector.hpp"

#include "topology/TopologicalDOF.hpp"

namespace residuum {
	namespace fem {
		namespace quantity {

			template<typename BackendT>
			class QuantityEvaluator {
			public:

				// reduces forms over the whole mesh domain
				template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename QuantityFormsT, typename QuadratureT>
				requires evaluator::EvalModel<ModelT, EvalQPT>
				static void evaluateDomain(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const QuantityFormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, BackendT>& U, const std::array<const linalg::types::Vector<Real, BackendT>*, EvalQPT::NumAuxStates>& auxStates, Real* out);

				// reduces forms over combo boundary faces (pointwise/strong flux)
				template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointBoundary EvalQPT, typename ModelT, typename QuantityFormsT, typename QuadratureT>
				requires evaluator::EvalModel<ModelT, EvalQPT>
				static void evaluateBoundaryRegistry(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const QuantityFormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, BackendT>& U, BoundaryQuantityRegistry<QuantityFormsT::TotalComponents>& registry);

				// reduces forms over combo boundary faces (consistent/weak flux): volume quadrature
				// accumulates the element vector, boundary quadrature accumulates measureSum. Called
				// once per distinct gather source an equation's residual needs (e.g. stiffness+source
				// with T, then mass with Tdot) -- zeroRegistry/computeMeasure should be true only for
				// the first such call, so later calls accumulate without re-zeroing or double-counting
				// the boundary measure (see the state-vector-unification issue: this two-pass-per-extra-
				// gather-source shape is a stopgap, not the intended long-term design)
				template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPVolT, evaluator::EvalQuadraturePointBoundary EvalQPBdyT, typename ModelT, typename FormsT, typename QuadratureVolumeT, typename QuadratureBoundaryT>
				requires evaluator::EvalModel<ModelT, EvalQPVolT>
				static void evaluateBoundaryReaction(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const FormsT& forms, const EvalEleT& evalEle, const QuadratureVolumeT& quadratureVolume, const QuadratureBoundaryT& quadratureBoundary, const linalg::types::Vector<Real, BackendT>& gatherSource, const std::array<const linalg::types::Vector<Real, BackendT>*, EvalQPVolT::NumAuxStates>& auxStates, bool zeroRegistry, bool computeMeasure, BoundaryQuantityRegistry<numDOFs>& registry);

				// recovers FieldT (e.g. a derived quantity like heat flux, or a pointwise material
				// property) to every mesh node, via RecoveryT (see fem/quantity/recovery); nodalOut has
				// length mesh.numNodes * FieldT::NumComponents
				template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FieldT, typename QuadratureT, typename RecoveryT>
				requires evaluator::EvalModel<ModelT, EvalQPT> && recovery::NodalRecoveryStrategy<RecoveryT>
				static void evaluateNodalField(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const FieldT& field, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, BackendT>& U, const RecoveryT& recovery, Real* nodalOut);

			}; // class QuantityEvaluator

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#include "backend/cpu/QuantityEvaluator.tpp"

#endif
