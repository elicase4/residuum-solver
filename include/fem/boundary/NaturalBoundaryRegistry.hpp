#ifndef RESIDUUM_FEM_BOUNDARY_NATURALBOUNDARYREGISTRY_HPP
#define RESIDUUM_FEM_BOUNDARY_NATURALBOUNDARYREGISTRY_HPP

#include <memory>
#include <unordered_map>
#include <vector>

#include "BoundaryCondition.hpp"
#include "NaturalBoundaryOperator.hpp"

#include "fem/evaluator/EvalModel.hpp"
#include "fem/evaluator/EvalQuadraturePointBoundary.hpp"

namespace residuum {
	namespace fem {
		namespace boundary {

			template<evaluator::EvalQuadraturePointBoundary EvalQPT>
			class NaturalBoundaryRegistry {

			public:

				template<BoundaryFunction FunctionT, typename FormsT, typename ModelT>
				requires evaluator::EvalModel<ModelT, EvalQPT>
				void registerBC(std::shared_ptr<BoundaryCondition<FunctionT>> bc, FormsT& forms, ModelT& model) {
					entries_[bc->tag].push_back(std::make_unique<NaturalBoundaryOperator<EvalQPT, FunctionT, FormsT, ModelT>>(std::move(bc), forms, model));
				}

				void apply(Int tag, EvalQPT& qp, Real* Fe) const {

					const auto* entries = getEntries(tag);

					if (!entries)
						return;

					for (const auto& entry : *entries)
						entry->apply(qp, Fe);
				}

				void gatherFaceElementData(Int tag, const Index* faceNodeGlobalIDs, Index nodesPerFace) const {

					const auto* entries = getEntries(tag);

					if (!entries)
						return;

					for (const auto& entry : *entries)
						entry->gatherFaceElementData(faceNodeGlobalIDs, nodesPerFace);
				}

				bool hasAny(Int tag) const {
					return entries_.contains(tag);
				}

				const std::vector<std::unique_ptr<NaturalBoundaryOperatorBase<EvalQPT>>>* getEntries(Int tag) const {

					auto it = entries_.find(tag);

					if (it == entries_.end())
						return nullptr;

					return &it->second;
				}
				
				bool isNatural(Int tag, Index component) const {
					
					const auto* entries = getEntries(tag);
					
					if (!entries) {
						return false;
					}

					for (const auto& entry : *entries) {
						if (entry->componentType(component) == BCCategory::Natural){
							return true;
						}
					}

					return false;
				
				}


			private:

				std::unordered_map<Int, std::vector<std::unique_ptr<NaturalBoundaryOperatorBase<EvalQPT>>>> entries_;

			}; // class NaturalBoundaryRegistry

		} // namespace boundary
	} // namespace fem
} // namespace residuum

#endif
