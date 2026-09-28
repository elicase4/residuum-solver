#ifndef RESIDUUM_FEM_QUADRATURE_GAUSSQUADRATUREHEX_HPP
#define RESIDUUM_FEM_QUADRATURE_GAUSSQUADRATUREHEX_HPP

#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/quadrature/GaussQuadrature1D.hpp"

namespace residuum {
	namespace fem {
		namespace quadrature {

			class GaussQuadratureHex {
			public:

				GaussQuadratureHex(Index numPointsX, Index numPointsY, Index numPointsZ) : quadX_(numPointsX), quadY_(numPointsY), quadZ_(numPointsZ), numPointsTotal_(numPointsX * numPointsY * numPointsZ) {}

				PDE_HOST PDE_DEVICE PDE_INLINE void getPoints(Real* xi) const;
				PDE_HOST PDE_DEVICE PDE_INLINE void getWeights(Real* w) const;

				Index numPointsXi() const { return quadX_.numPoints(); }
				Index numPointsEta() const { return quadY_.numPoints(); }
				Index numPointsZeta() const { return quadZ_.numPoints(); }
				Index numPointsTotal() const { return numPointsTotal_; }

			private:
				GaussQuadrature1D quadX_;
				GaussQuadrature1D quadY_;
				GaussQuadrature1D quadZ_;
				Index numPointsTotal_;

			}; // class GaussQuadratureHex
			
		} // namespace quadrature
	} // namespace fem
} // namespace residuum

#include "GaussQuadratureHex.tpp"

#endif
