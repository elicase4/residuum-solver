#ifndef RESIDUUM_FEM_QUANTITY_BOUNDARYQUADRATUREMONITORGROUP_HPP
#define RESIDUUM_FEM_QUANTITY_BOUNDARYQUADRATUREMONITORGROUP_HPP

#include <utility>

#include "core/Types.hpp"

#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/evaluator/EvalElement.hpp"
#include "fem/evaluator/EvalModel.hpp"
#include "fem/evaluator/EvalQuadraturePointBoundary.hpp"
#include "fem/quantity/MonitorAccumulator.hpp"
#include "fem/quantity/MonitorGroup.hpp"
#include "fem/quantity/QuantityEvaluator.hpp"

#include "linalg/types/Vector.hpp"

#include "mesh/Mesh.hpp"

#include "topology/TopologicalDOF.hpp"

namespace residuum {
	namespace fem {
		namespace quantity {

			// pointwise (strong) boundary flux: a QuantityForm evaluated at boundary quadrature points
			template<typename BackendT, Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointBoundary EvalQPT, typename ModelT, typename QuantityFormsT, typename QuadratureT>
			requires evaluator::EvalModel<ModelT, EvalQPT>
			class BoundaryQuadratureMonitorGroup : public MonitorGroup {
			public:

				BoundaryQuadratureMonitorGroup(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const ModelT& model, QuantityFormsT forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, BackendT>& U)
					: mesh_(mesh), topoDOF_(topoDOF), bcRegistry_(bcRegistry), model_(model), forms_(std::move(forms)), evalEle_(evalEle), quadrature_(quadrature), U_(U) {}

				void registerTag(Int tag) override { accumulator_.registerTag(tag); }

				Index addCombination() override { return accumulator_.addCombination(); }

				void addTerm(Index combinationIndex, Int tag, Real coefficient) override { accumulator_.addTerm(combinationIndex, tag, coefficient); }

				void evaluate(Real time, const std::function<void(Index, const Real*)>& callback) override {
					QuantityEvaluator<BackendT>::template evaluateBoundaryRegistry<numDOFs, EvalEleT, EvalQPT, ModelT, QuantityFormsT, QuadratureT>(mesh_, topoDOF_, bcRegistry_, time, model_, forms_, evalEle_, quadrature_, U_, accumulator_.registry());
					accumulator_.forEachCombination(callback);
				}

			private:

				const mesh::Mesh& mesh_;
				const topology::TopologicalDOF<numDOFs>& topoDOF_;
				const fem::boundary::EssentialBoundaryRegistry& bcRegistry_;
				const ModelT& model_;
				QuantityFormsT forms_;
				const EvalEleT& evalEle_;
				const QuadratureT& quadrature_;
				const linalg::types::Vector<Real, BackendT>& U_;

				MonitorAccumulator<QuantityFormsT::TotalComponents> accumulator_;

			}; // class BoundaryQuadratureMonitorGroup

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
