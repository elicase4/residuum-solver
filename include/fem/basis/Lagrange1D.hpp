#ifndef RESIDUUM_FEM_BASIS_LAGRANGE1D_HPP
#define RESIDUUM_FEM_BASIS_LAGRANGE1D_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace residuum {
	namespace fem {
		namespace basis {

			class Lagrange1D {
			public:

				explicit Lagrange1D(Index order) : order_(order) {}

				PDE_HOST PDE_DEVICE PDE_INLINE void eval(Real xi, Real* N) const;
				PDE_HOST PDE_DEVICE PDE_INLINE void evalFirstDerivative(Real xi, Real* dN) const;
				PDE_HOST PDE_DEVICE PDE_INLINE void evalSecondDerivative(Real xi, Real* d2N) const;

				Index order() const { return order_; }
				Index nodesPerElement() const { return order_ + 1; }

				static constexpr Index ParametricDim = 1;

			private:
				Index order_;

			}; // class Lagrange1D

		} // namespace basis
	} // namespace fem
} // namespace residuum

#include "Lagrange1D.tpp"

#endif
