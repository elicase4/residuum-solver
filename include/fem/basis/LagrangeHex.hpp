#ifndef RESIDUUM_FEM_BASIS_LAGRANGEHEX_HPP
#define RESIDUUM_FEM_BASIS_LAGRANGEHEX_HPP

#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/basis/Lagrange1D.hpp"

namespace residuum {
	namespace fem {
		namespace basis {

			class LagrangeHex {
			public:

				LagrangeHex(Index px, Index py, Index pz) : basisX_(px), basisY_(py), basisZ_(pz), nodesPerElement_((px + 1) * (py + 1) * (pz + 1)) {}

				PDE_HOST PDE_DEVICE PDE_INLINE void eval(const Real* xi, Real* N) const;
				PDE_HOST PDE_DEVICE PDE_INLINE void evalGradient(const Real* xi, Real* dNdxi) const;
				PDE_HOST PDE_DEVICE PDE_INLINE void evalHessian(const Real* xi, Real* d2Nd2xi) const;
				PDE_HOST PDE_DEVICE PDE_INLINE void evalLaplacian(const Real* xi, Real* lapN) const;

				PDE_HOST PDE_DEVICE PDE_INLINE void getFaceTopology(const Int rngID, Real* nRef) const;
				PDE_HOST PDE_DEVICE PDE_INLINE Index nodesPerFace(const Int rngID) const;
				PDE_HOST PDE_DEVICE PDE_INLINE void getFaceNodes(const Int rngID, Index* nodeIDs) const;
				PDE_HOST PDE_DEVICE PDE_INLINE void mapFaceToElement(const Int rngID, const Real* xi_face, Real* xi_elem) const;

				Index orderX() const { return basisX_.order(); }
				Index orderY() const { return basisY_.order(); }
				Index orderZ() const { return basisZ_.order(); }
				Index nodesPerElement() const { return nodesPerElement_; }

				static constexpr Index ParametricDim = 3;

			private:
				Lagrange1D basisX_;
				Lagrange1D basisY_;
				Lagrange1D basisZ_;
				Index nodesPerElement_;

			}; // class LagrangeHex

		} // namespace basis
	} // namespace fem
} // namespace residuum

#include "LagrangeHex.tpp"

#endif
