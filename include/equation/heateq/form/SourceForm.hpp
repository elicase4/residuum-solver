#ifndef RESIDUUM_EQUATION_HEATEQ_FORM_SOURCEFORM_HPP
#define RESIDUUM_EQUATION_HEATEQ_FORM_SOURCEFORM_HPP

#include <type_traits>
#include <utility>

#include "fem/form/LinearForm.hpp"
#include "equation/heateq/evaluator/SourceFunction.hpp"

namespace residuum::equation::heateq {

	template<typename QuadraturePointVolumeT, typename SourceFunctionT>
	requires fem::evaluator::EvalFunction<SourceFunctionT>
	struct SourceForm {
		
		SourceFunctionT sourceFunction;
		
		template<typename... Args>
		requires (!(sizeof...(Args) == 1 && (std::is_same_v<std::remove_cvref_t<Args>, SourceForm> && ...)))
		constexpr SourceForm(Args&&... args) : sourceFunction(std::forward<Args>(args)...) {}

		PDE_HOST PDE_DEVICE void computeElementLevelVector(const QuadraturePointVolumeT& qp, const Real*, Real* Fe) const {
			
			Real val[SourceFunctionT::NumComponents];
			sourceFunction.eval(qp.time, qp.x, val);

			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				for (Index i = 0; i < SourceFunctionT::NumComponents; ++i) {
					Fe[a + i] += (val[i] * qp.N[a]) * qp.measure * qp.w;
				}
			}
		}

	};// struct SourceForm

} // namespace residuum::equation::heateq

#endif
