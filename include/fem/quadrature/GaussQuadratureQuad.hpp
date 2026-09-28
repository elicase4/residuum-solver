#ifndef RESIDUUM_FEM_QUADRATURE_GAUSSQUADRATUREQUAD_HPP
#define RESIDUUM_FEM_QUADRATURE_GAUSSQUADRATUREQUAD_HPP

#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/quadrature/GaussQuadrature1D.hpp"

namespace residuum {
	namespace fem {
		namespace quadrature {

			class GaussQuadratureQuad {
			public:

				GaussQuadratureQuad(Index numPointsX, Index numPointsY) : quadX_(numPointsX), quadY_(numPointsY), numPointsTotal_(numPointsX * numPointsY) {}

				PDE_HOST PDE_DEVICE PDE_INLINE void getPoints(Real* xi) const;
				PDE_HOST PDE_DEVICE PDE_INLINE void getWeights(Real* w) const;

				Index numPointsXi() const { return quadX_.numPoints(); }
				Index numPointsEta() const { return quadY_.numPoints(); }
				Index numPointsTotal() const { return numPointsTotal_; }

			private:
				GaussQuadrature1D quadX_;
				GaussQuadrature1D quadY_;
				Index numPointsTotal_;

			}; // class GaussQuadratureQuad
			
		} // namespace quadrature
	} // namespace fem
} // namespace residuum

#include "GaussQuadratureQuad.tpp"

#endif
