#ifndef RESIDUUM_FEM_QUANTITY_QUANTITYEVALUATOR_HPP
#define RESIDUUM_FEM_QUANTITY_QUANTITYEVALUATOR_HPP

#include "config/Platform.hpp"

#include "core/Types.hpp"

#include "fem/assembly/Assembler.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/evaluator/EvalElement.hpp"
#include "fem/evaluator/EvalModel.hpp"
#include "fem/evaluator/EvalQuadraturePointBoundary.hpp"
#include "fem/evaluator/EvalQuadraturePointVolume.hpp"
#include "fem/evaluator/StateVector.hpp"
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

				template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointBoundary EvalQPT, typename ModelT, typename QuantityFormsT, typename QuadratureT, typename StateVectorT>
				requires evaluator::EvalModel<ModelT, EvalQPT>
				static void evaluateBoundaryRegistry(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const QuantityFormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, StateVectorT& state, BoundaryQuantityRegistry<QuantityFormsT::TotalComponents>& registry);

				template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPVolT, evaluator::EvalQuadraturePointBoundary EvalQPBdyT, typename ModelT, typename QuadratureVolumeT, typename QuadratureBoundaryT, typename StateVectorT, typename... FormsWithOperandT>
				requires evaluator::EvalModel<ModelT, EvalQPVolT>
				static void evaluateBoundaryReaction(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const EvalEleT& evalEle, const QuadratureVolumeT& quadratureVolume, const QuadratureBoundaryT& quadratureBoundary, StateVectorT& state, BoundaryQuantityRegistry<numDOFs>& registry, FormsWithOperandT... pairs);

				template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FieldT, typename QuadratureT, typename RecoveryT, typename StateVectorT>
				requires evaluator::EvalModel<ModelT, EvalQPT> && recovery::NodalRecoveryStrategy<RecoveryT>
				static void evaluateNodalField(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const FieldT& field, const EvalEleT& evalEle, const QuadratureT& quadrature, StateVectorT& state, const RecoveryT& recovery, Real* nodalOut);

			}; // class QuantityEvaluator

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#include "backend/cpu/QuantityEvaluator.tpp"

#endif
