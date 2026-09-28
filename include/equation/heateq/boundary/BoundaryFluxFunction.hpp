#ifndef RESIDUUM_EQUATION_HEATEQ_BOUNDARY_BOUNDARYFLUXFUNCTION_HPP
#define RESIDUUM_EQUATION_HEATEQ_BOUNDARY_BOUNDARYFLUXFUNCTION_HPP

#include <cmath>
#include <type_traits>
#include <utility>

#include "core/Types.hpp"
#include "fem/evaluator/EvalFunction.hpp"

namespace residuum::equation::heateq {

	template<Index SpatialDimension, Index numDOFs, class CallableT>
	struct BoundaryFluxFunction {

		static constexpr Index NumComponents = numDOFs;
		static constexpr Index SpatialDim = SpatialDimension;

		CallableT f;

		template<typename... Args>
		requires (!(sizeof...(Args) == 1 && (std::is_same_v<std::remove_cvref_t<Args>, BoundaryFluxFunction> && ...)))
		constexpr BoundaryFluxFunction(Args&&... args) : f(std::forward<Args>(args)...) {}

		void eval(const Real time, const Real* x, Real* outValue) const {
			f(time, x, outValue);
		}

		void evalGradient(const Real, const Real*, Real*) const {}

	}; // struct BoundaryFluxFunction

} // namespace residuum::equation::heateq

#endif
