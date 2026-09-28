#ifndef RESIDUUM_EQUATION_HEATEQ_FORM_FLUXBOUNDARYFORM_HPP
#define RESIDUUM_EQUATION_HEATEQ_FORM_FLUXBOUNDARYFORM_HPP

#include <type_traits>
#include <utility>

#include "fem/boundary/BoundaryCondition.hpp"
#include "fem/evaluator/EvalQuadraturePointBoundary.hpp"
#include "fem/form/LinearForm.hpp"

namespace residuum::equation::heateq {

	template<typename QuadraturePointBoundaryT, typename FluxFunctionT>
	struct FluxBoundaryForm {

		FluxFunctionT fluxFunction;
		
		template<typename... Args>
		requires (!(sizeof...(Args) == 1 && (std::is_same_v<std::remove_cvref_t<Args>, FluxBoundaryForm> && ...)))
		constexpr FluxBoundaryForm(Args&&... args) : fluxFunction(std::forward<Args>(args)...) {}
		
		PDE_HOST PDE_DEVICE void computeElementLevelVector(const QuadraturePointBoundaryT& qp, const Real*, Real* Fe) const {
		
			Real val[FluxFunctionT::SpatialDim * FluxFunctionT::NumComponents];
			fluxFunction.eval(qp.time, qp.x, val);

			for (Index a = 0; a < qp.nodesPerFace(); ++a){
				for (Index i = 0; i < FluxFunctionT::NumComponents; ++i) {
					for (Index sD = 0; sD < FluxFunctionT::SpatialDim; ++sD){
						Fe[a*FluxFunctionT::NumComponents + i] += (val[i*FluxFunctionT::SpatialDim + sD] * qp.Nface[a]) * qp.normal[sD] * qp.w;
					}
				}
			}

		}

	};// struct FluxBoundaryForm

} // namespace residuum::equation::heateq

#endif
