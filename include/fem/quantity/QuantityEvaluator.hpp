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

				// reduces forms over combo boundary faces
				template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointBoundary EvalQPT, typename ModelT, typename QuantityFormsT, typename QuadratureT>
				requires evaluator::EvalModel<ModelT, EvalQPT>
				static void evaluateBoundaryRegistry(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const QuantityFormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, BackendT>& U, BoundaryQuantityRegistry<QuantityFormsT>& registry);

			}; // class QuantityEvaluator

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#include "backend/cpu/QuantityEvaluator.tpp"

#endif
