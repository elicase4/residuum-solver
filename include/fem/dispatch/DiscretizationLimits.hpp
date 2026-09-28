#ifndef RESIDUUM_FEM_DISPATCH_DISCRETIZATIONLIMITS_HPP
#define RESIDUUM_FEM_DISPATCH_DISCRETIZATIONLIMITS_HPP

#include "core/Types.hpp"

namespace residuum {
	namespace fem {
		namespace dispatch {

			inline constexpr Index kMaxBasisOrder = 3;
			inline constexpr Index kMaxQuadraturePoints1D = 5;

			namespace detail {

				constexpr Index ipow(Index base, Index exp) {
					Index result = 1;
					for (Index i = 0; i < exp; ++i) {
						result *= base;
					}
					return result;
				}

			} // namespace detail

			template<Index NPD>
			inline constexpr Index kMaxNodesPerElement = detail::ipow(kMaxBasisOrder + 1, NPD);

			template<Index NPD>
			inline constexpr Index kMaxNodesPerElementBoundary = detail::ipow(kMaxBasisOrder + 1, NPD-1);

			template<Index NPD>
			inline constexpr Index kMaxQuadraturePointsTotal = detail::ipow(kMaxQuadraturePoints1D, NPD);

			template<Index NPD>
			inline constexpr Index kMaxQuadraturePointsTotalBoundary = detail::ipow(kMaxQuadraturePoints1D, NPD - 1);

		} // namespace dispatch
	} // namespace fem
} // namespace residuum

#endif
