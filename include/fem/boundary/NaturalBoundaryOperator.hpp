#ifndef RESIDUUM_FEM_BOUNDARY_NATURALBOUNDARYOPERATOR_HPP
#define RESIDUUM_FEM_BOUNDARY_NATURALBOUNDARYOPERATOR_HPP

#include <memory>
#include <unordered_map>
#include <vector>

#include "BoundaryCondition.hpp"

#include "fem/evaluator/EvalModel.hpp"
#include "fem/evaluator/EvalQuadraturePointBoundary.hpp"

namespace residuum {
	namespace fem {
		namespace boundary {

			template<evaluator::EvalQuadraturePointBoundary EvalQPT>
			class NaturalBoundaryOperatorBase {
			public:

				virtual ~NaturalBoundaryOperatorBase() = default;
				virtual void apply(EvalQPT&, Real*) const = 0;
				virtual BCCategory componentType(Index component) const = 0;

				virtual void gatherFaceElementData(const Index*, Index) const {}

			}; // class NaturalBoundaryOperatorBase

			template<evaluator::EvalQuadraturePointBoundary EvalQPT, BoundaryFunction FunctionT, typename FormsT, typename ModelT>
			requires evaluator::EvalModel<ModelT, EvalQPT>
			class NaturalBoundaryOperator final : public NaturalBoundaryOperatorBase<EvalQPT> {
			public:

				NaturalBoundaryOperator(std::shared_ptr<BoundaryCondition<FunctionT>> bc, FormsT& forms, ModelT& model) : bc_(std::move(bc)), forms_(forms), model_(model) {}

				void apply(EvalQPT& qp, Real* Fe) const override {
					model_.eval(qp);
					model_.evalGradient(qp);
					forms_.computeElementLevelVector(qp, nullptr, Fe);
				}

				BCCategory componentType(Index component) const override {
					return bc_->componentType[component];
				}

				void gatherFaceElementData(const Index* faceNodeGlobalIDs, Index nodesPerFace) const override {
					forms_.gatherElementData(faceNodeGlobalIDs, nodesPerFace);
				}

			private:

				std::shared_ptr<BoundaryCondition<FunctionT>> bc_;
				FormsT& forms_;
				ModelT& model_;

			}; // class NaturalBoundaryOperator

		} // namespace boundary
	} // namespace fem
} // namespace residuum

#endif
