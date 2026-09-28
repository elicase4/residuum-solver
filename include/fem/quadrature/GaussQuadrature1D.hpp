#ifndef RESIDUUM_FEM_QUADRATURE_GAUSSQUADRATURE1D_HPP
#define RESIDUUM_FEM_QUADRATURE_GAUSSQUADRATURE1D_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace residuum {
	namespace fem {
		namespace quadrature {

			class GaussQuadrature1D {
			public:

				explicit GaussQuadrature1D(Index numPoints) : numPoints_(numPoints) {}

				PDE_HOST PDE_DEVICE PDE_INLINE void getPoints(Real* xi) const;
				PDE_HOST PDE_DEVICE PDE_INLINE void getWeights(Real* w) const;

				Index numPoints() const { return numPoints_; }
				Index numPointsTotal() const { return numPoints_; }

			private:
				Index numPoints_;

			}; // class GaussQuadrature1D

		} // namespace quadrature
	} // namespace fem
} // namespace residuum

#include "GaussQuadrature1D.tpp"

#endif
