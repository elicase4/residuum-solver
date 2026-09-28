#ifndef RESIDUUM_FEM_QUANTITY_MONITORGROUPIMPL_HPP
#define RESIDUUM_FEM_QUANTITY_MONITORGROUPIMPL_HPP

#include <utility>
#include <vector>

#include "core/Types.hpp"

#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/evaluator/EvalElement.hpp"
#include "fem/evaluator/EvalModel.hpp"
#include "fem/evaluator/EvalQuadraturePointBoundary.hpp"
#include "fem/quantity/BoundaryQuantityCombination.hpp"
#include "fem/quantity/BoundaryQuantityRegistry.hpp"
#include "fem/quantity/MonitorGroup.hpp"
#include "fem/quantity/QuantityEvaluator.hpp"

#include "linalg/types/Vector.hpp"

#include "mesh/Mesh.hpp"

#include "topology/TopologicalDOF.hpp"

namespace residuum {
	namespace fem {
		namespace quantity {

			template<typename BackendT, Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointBoundary EvalQPT, typename ModelT, typename QuantityFormsT, typename QuadratureT>
			requires evaluator::EvalModel<ModelT, EvalQPT>
			class MonitorGroupImpl : public MonitorGroup {
			public:

				MonitorGroupImpl(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const ModelT& model, QuantityFormsT forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, BackendT>& U) : mesh_(mesh), topoDOF_(topoDOF), bcRegistry_(bcRegistry), model_(model), forms_(std::move(forms)), evalEle_(evalEle), quadrature_(quadrature), U_(U) {}

				void registerTag(Int tag) override { registry_.registerTag(tag); }

				Index addCombination() override { combinations_.emplace_back(); return static_cast<Index>(combinations_.size() - 1); }

				void addTerm(Index combinationIndex, Int tag, Real coefficient) override { combinations_[combinationIndex].addTerm(tag, coefficient); }

				void evaluate(Real time, const std::function<void(Index, const Real*)>& sink) override {

					QuantityEvaluator<BackendT>::template evaluateBoundaryRegistry<numDOFs, EvalEleT, EvalQPT, ModelT, QuantityFormsT, QuadratureT>(mesh_, topoDOF_, bcRegistry_, time, model_, forms_, evalEle_, quadrature_, U_, registry_);

					Real value[QuantityFormsT::TotalComponents];
					for (Index i = 0; i < static_cast<Index>(combinations_.size()); ++i) {
						combinations_[i].evaluate(registry_, value);
						sink(i, value);
					}

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

				BoundaryQuantityRegistry<QuantityFormsT> registry_;
				std::vector<BoundaryQuantityCombination<QuantityFormsT>> combinations_;

			}; // class MonitorGroupImpl

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
