#ifndef RESIDUUM_FEM_QUANTITY_BOUNDARYREACTIONMONITORGROUP_HPP
#define RESIDUUM_FEM_QUANTITY_BOUNDARYREACTIONMONITORGROUP_HPP

#include <algorithm>
#include <array>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "core/Types.hpp"

#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/evaluator/EvalElement.hpp"
#include "fem/evaluator/EvalModel.hpp"
#include "fem/evaluator/EvalQuadraturePointBoundary.hpp"
#include "fem/evaluator/EvalQuadraturePointVolume.hpp"
#include "fem/quantity/BoundaryVectorProjection.hpp"
#include "fem/quantity/MonitorAccumulator.hpp"
#include "fem/quantity/MonitorGroup.hpp"
#include "fem/quantity/QuantityEvaluator.hpp"
#include "fem/quantity/Reduction.hpp"

#include "linalg/types/Vector.hpp"

#include "mesh/Mesh.hpp"

#include "topology/TopologicalDOF.hpp"

namespace residuum {
	namespace fem {
		namespace quantity {

			// weak/consistent flux: FormsT (a FormRegistry bundling whatever residual contributions the
			// equation needs with a shared gather source, e.g. stiffness+source with T) is volume-
			// assembled per element touching a registered tag, reduced to that tag's face-local rows;
			// boundary quadrature over the same faces tracks measureSum (see
			// QuantityEvaluator::evaluateBoundaryReaction). An optional second contribution (FormsT2,
			// e.g. mass with Tdot) accumulates on top of the first for equations with more than one
			// distinct gather source -- a stopgap until the state-vector-unification issue lands, since
			// it costs a second full element-touching-a-tag pass rather than a single fused one
			template<typename BackendT, Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPVolT, evaluator::EvalQuadraturePointBoundary EvalQPBdyT, typename ModelT, typename FormsT, typename QuadratureVolumeT, typename QuadratureBoundaryT, typename FormsT2 = FormsT, typename ModelT2 = ModelT>
			requires evaluator::EvalModel<ModelT, EvalQPVolT>
			class BoundaryReactionMonitorGroup : public MonitorGroup {
			public:

				struct SecondContribution {
					FormsT2 forms;
					const ModelT2& model;
					const linalg::types::Vector<Real, BackendT>& gatherSource;
					std::array<const linalg::types::Vector<Real, BackendT>*, EvalQPVolT::NumAuxStates> auxStates;
				}; // struct SecondContribution

				BoundaryReactionMonitorGroup(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const ModelT& model, FormsT forms, const EvalEleT& evalEle, const QuadratureVolumeT& quadratureVolume, const QuadratureBoundaryT& quadratureBoundary, const linalg::types::Vector<Real, BackendT>& U, const std::array<const linalg::types::Vector<Real, BackendT>*, EvalQPVolT::NumAuxStates>& auxStates, Reduction reduction, std::optional<SecondContribution> second = std::nullopt)
					: mesh_(mesh), topoDOF_(topoDOF), bcRegistry_(bcRegistry), model_(model), forms_(std::move(forms)), evalEle_(evalEle), quadratureVolume_(quadratureVolume), quadratureBoundary_(quadratureBoundary), U_(U), auxStates_(auxStates), reduction_(reduction), second_(std::move(second)) {}

				void registerTag(Int tag) override { accumulator_.registerTag(tag); }

				Index addCombination() override { return accumulator_.addCombination(); }

				void addTerm(Index combinationIndex, Int tag, Real coefficient) override { accumulator_.addTerm(combinationIndex, tag, coefficient); }

				// separate index space from addCombination/addTerm; e.g. lift as a projection of the
				// force vector this monitor already produces, rather than its raw components
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

					QuantityEvaluator<BackendT>::template evaluateBoundaryReaction<numDOFs, EvalEleT, EvalQPVolT, EvalQPBdyT, ModelT, FormsT, QuadratureVolumeT, QuadratureBoundaryT>(mesh_, topoDOF_, bcRegistry_, time, model_, forms_, evalEle_, quadratureVolume_, quadratureBoundary_, U_, auxStates_, /*zeroRegistry=*/true, /*computeMeasure=*/true, accumulator_.registry());

					if (second_) {
						QuantityEvaluator<BackendT>::template evaluateBoundaryReaction<numDOFs, EvalEleT, EvalQPVolT, EvalQPBdyT, ModelT2, FormsT2, QuadratureVolumeT, QuadratureBoundaryT>(mesh_, topoDOF_, bcRegistry_, time, second_->model, second_->forms, evalEle_, quadratureVolume_, quadratureBoundary_, second_->gatherSource, second_->auxStates, /*zeroRegistry=*/false, /*computeMeasure=*/false, accumulator_.registry());
					}

					// FormsT has no Reduction concept of its own (unlike QuantityForms), so this monitor
					// applies Average-mode division itself, on the one quantity it ever produces
					if (reduction_ == Reduction::Average) {
						for (auto& [tag, entry] : accumulator_.registry().entries()) {
							for (Index c = 0; c < numDOFs; ++c) entry.result[c] /= entry.measureSum;
						}
					}

					// a monitor is configured as either combinations or weighted projections, never
					// both, so exactly one of these two loops ever has anything to report
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
				FormsT forms_;
				const EvalEleT& evalEle_;
				const QuadratureVolumeT& quadratureVolume_;
				const QuadratureBoundaryT& quadratureBoundary_;
				const linalg::types::Vector<Real, BackendT>& U_;
				std::array<const linalg::types::Vector<Real, BackendT>*, EvalQPVolT::NumAuxStates> auxStates_;
				Reduction reduction_;
				std::optional<SecondContribution> second_;

				MonitorAccumulator<numDOFs> accumulator_;
				std::vector<BoundaryVectorProjection<numDOFs>> projections_;

			}; // class BoundaryReactionMonitorGroup

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
