#ifndef RESIDUUM_FEM_QUANTITY_BOUNDARYREACTIONMONITORGROUP_HPP
#define RESIDUUM_FEM_QUANTITY_BOUNDARYREACTIONMONITORGROUP_HPP

#include <array>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "core/Types.hpp"

#include "fem/assembly/Assembler.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/evaluator/EvalElement.hpp"
#include "fem/evaluator/EvalModel.hpp"
#include "fem/evaluator/EvalQuadraturePointBoundary.hpp"
#include "fem/evaluator/EvalQuadraturePointVolume.hpp"
#include "fem/evaluator/StateVector.hpp"
#include "fem/quantity/BoundaryVectorProjection.hpp"
#include "fem/quantity/MonitorAccumulator.hpp"
#include "fem/quantity/MonitorGroup.hpp"
#include "fem/quantity/QuantityEvaluator.hpp"
#include "fem/quantity/Reduction.hpp"

#include "mesh/Mesh.hpp"

#include "topology/TopologicalDOF.hpp"

namespace residuum {
	namespace fem {
		namespace quantity {

			template<typename BackendT, Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPVolT, evaluator::EvalQuadraturePointBoundary EvalQPBdyT, typename ModelT, typename QuadratureVolumeT, typename QuadratureBoundaryT, typename StateVectorT, typename... FormsWithOperandT>
			requires evaluator::EvalModel<ModelT, EvalQPVolT>
			class BoundaryReactionMonitorGroup : public MonitorGroup {
			public:

				BoundaryReactionMonitorGroup(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const ModelT& model, const EvalEleT& evalEle, const QuadratureVolumeT& quadratureVolume, const QuadratureBoundaryT& quadratureBoundary, StateVectorT& state, Reduction reduction, FormsWithOperandT... pairs)
					: mesh_(mesh), topoDOF_(topoDOF), bcRegistry_(bcRegistry), model_(model), evalEle_(evalEle), quadratureVolume_(quadratureVolume), quadratureBoundary_(quadratureBoundary), state_(state), reduction_(reduction), pairs_(pairs...) {}

				void registerTag(Int tag) override { accumulator_.registerTag(tag); }

				Index addCombination() override { return accumulator_.addCombination(); }

				void addTerm(Index combinationIndex, Int tag, Real coefficient) override { accumulator_.addTerm(combinationIndex, tag, coefficient); }

				Index addWeightedCombination() override { projections_.emplace_back(); return static_cast<Index>(projections_.size() - 1); }

				void addWeightedTerm(Index combinationIndex, Int tag, const std::vector<Real>& weights) override {
					if (weights.size() != numDOFs) {
						throw std::runtime_error("BoundaryReactionMonitorGroup: weights.size() (" + std::to_string(weights.size()) + ") must equal numDOFs (" + std::to_string(numDOFs) + ")");
					}
					std::array<Real, numDOFs> w;
					std::copy(weights.begin(), weights.end(), w.begin());
					projections_[combinationIndex].addTerm(tag, w);
				}

				void evaluate(Real time, const std::function<void(Index, const Real*)>& callback) override {

					std::apply([&](auto&... pairs){
						QuantityEvaluator<BackendT>::template evaluateBoundaryReaction<numDOFs, EvalEleT, EvalQPVolT, EvalQPBdyT, ModelT, QuadratureVolumeT, QuadratureBoundaryT, StateVectorT, FormsWithOperandT...>(mesh_, topoDOF_, bcRegistry_, time, model_, evalEle_, quadratureVolume_, quadratureBoundary_, state_, accumulator_.registry(), pairs...);
					}, pairs_);

					// apply Average-mode division
					if (reduction_ == Reduction::Average) {
						for (auto& [tag, entry] : accumulator_.registry().entries()) {
							for (Index c = 0; c < numDOFs; ++c) entry.result[c] /= entry.measureSum;
						}
					}

					// a monitor is configured as either combinations or weighted projections
					accumulator_.forEachCombination(callback);
					for (Index i = 0; i < static_cast<Index>(projections_.size()); ++i) {
						Real value = projections_[i].evaluate(accumulator_.registry());
						callback(i, &value);
					}

				}

			private:

				const mesh::Mesh& mesh_;
				const topology::TopologicalDOF<numDOFs>& topoDOF_;
				const fem::boundary::EssentialBoundaryRegistry& bcRegistry_;
				const ModelT& model_;
				const EvalEleT& evalEle_;
				const QuadratureVolumeT& quadratureVolume_;
				const QuadratureBoundaryT& quadratureBoundary_;
				StateVectorT& state_;
				Reduction reduction_;
				std::tuple<FormsWithOperandT...> pairs_;

				MonitorAccumulator<numDOFs> accumulator_;
				std::vector<BoundaryVectorProjection<numDOFs>> projections_;

			}; // class BoundaryReactionMonitorGroup

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
