#ifndef RESIDUUM_EQUATION_HEATEQ_FORM_NODALSOURCEFORM_HPP
#define RESIDUUM_EQUATION_HEATEQ_FORM_NODALSOURCEFORM_HPP

#include <utility>

#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/evaluator/EvalNodalData.hpp"
#include "fem/evaluator/EvalQuadraturePointVolume.hpp"
#include "fem/form/LinearForm.hpp"

namespace residuum::equation::heateq {

	template<typename QuadraturePointVolumeT, fem::evaluator::EvalNodalData SourceT, Index NumComponents_>
	class NodalSourceForm {
	public:

		static constexpr Index NumComponents = NumComponents_;

		template<typename... Args>
		requires (!(sizeof...(Args) == 1 && (std::is_same_v<std::remove_cvref_t<Args>, NodalSourceForm> && ...)))
		constexpr NodalSourceForm(Args&&... args) : source_(std::forward<Args>(args)...) {}

		// per-element, before the quadrature-point loop
		void gatherElementData(const Index* nodeIDs, Index nodesPerElement) const {
			for (Index a = 0; a < nodesPerElement; ++a) {
				source_.eval(nodeIDs[a], &De_[a * NumComponents]);
			}
		}

		PDE_HOST PDE_DEVICE void computeElementLevelVector(const QuadraturePointVolumeT& qp, const Real*, Real* Fe) const {

			Real val[NumComponents] = {0};

			// interpolate the gathered per-node source data to this QP
			for (Index b = 0; b < qp.nodesPerElement(); ++b) {
				for (Index i = 0; i < NumComponents; ++i) {
					val[i] += qp.N[b] * De_[b*NumComponents + i];
				}
			}

			// same weak-form accumulation as SourceForm
			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				for (Index i = 0; i < NumComponents; ++i) {
					Fe[a*NumComponents + i] += (val[i] * qp.N[a]) * qp.measure * qp.w;
				}
			}

		}

	private:

		SourceT source_;
		mutable Real De_[fem::dispatch::kMaxNodesPerElement<QuadraturePointVolumeT::ParametricDim> * NumComponents_];

	}; // class NodalSourceForm

} // namespace residuum::equation::heateq

#endif
