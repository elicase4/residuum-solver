#ifndef RESIDUUM_FEM_DISPATCH_DISCRETIZATIONDISPATCH_HPP
#define RESIDUUM_FEM_DISPATCH_DISCRETIZATIONDISPATCH_HPP

#include <stdexcept>
#include <string>
#include <utility>

#include "core/Types.hpp"

#include "fem/basis/LagrangeQuad.hpp"
#include "fem/basis/LagrangeHex.hpp"
#include "fem/dispatch/DiscretizationLimits.hpp"

#include "fem/quadrature/GaussQuadrature1D.hpp"
#include "fem/quadrature/GaussQuadratureQuad.hpp"
#include "fem/quadrature/GaussQuadratureHex.hpp"

#include "mesh/ElementFamily.hpp"

namespace residuum {
	namespace fem {
		namespace dispatch {

			template<mesh::ElementFamily Family> struct ElementTypeTraits;

			template<> struct ElementTypeTraits<mesh::ElementFamily::Quad> {
				using BasisType = basis::LagrangeQuad;
				using QuadratureVolumeType = quadrature::GaussQuadratureQuad;
				using QuadratureBoundaryType = quadrature::GaussQuadrature1D;
			};

			template<> struct ElementTypeTraits<mesh::ElementFamily::Hex> {
				using BasisType = basis::LagrangeHex;
				using QuadratureVolumeType = quadrature::GaussQuadratureHex;
				using QuadratureBoundaryType = quadrature::GaussQuadratureQuad;
			};

			// throws instead of silently letting an out-of-range order/quadrature-point count
			// reach the stack-allocated buffers in Assembler/BoundaryApplicator, which are sized
			// from these same limits and don't re-check them
			inline void validateDiscretizationLimits(Index basisOrder, Index quadraturePoints1D) {

				if (basisOrder < 1 || basisOrder > kMaxBasisOrder) {
					throw std::runtime_error("fem::dispatch: basis order " + std::to_string(basisOrder) + " is outside the supported range [1, " + std::to_string(kMaxBasisOrder) + "]");
				}

				if (quadraturePoints1D < 1 || quadraturePoints1D > kMaxQuadraturePoints1D) {
					throw std::runtime_error("fem::dispatch: quadrature point count " + std::to_string(quadraturePoints1D) + " is outside the supported range [1, " + std::to_string(kMaxQuadraturePoints1D) + "]");
				}

			}

			template<Index NSD, typename VisitorT>
			bool dispatchQuad(Index px, Index py, Index xi, Index eta, VisitorT&& visitor) {

				validateDiscretizationLimits(px, xi);
				validateDiscretizationLimits(py, eta);

				using Traits = ElementTypeTraits<mesh::ElementFamily::Quad>;

				typename Traits::BasisType basisInst(px, py);
				typename Traits::QuadratureVolumeType quadVolInst(xi, eta);

				const Index nBdy = (xi > eta) ? xi : eta;
				typename Traits::QuadratureBoundaryType quadBdyInst(nBdy);

				return visitor.template operator()<NSD, 2, mesh::ElementFamily::Quad>(basisInst, quadVolInst, quadBdyInst);

			}

			template<Index NSD, typename VisitorT>
			bool dispatchHex(Index px, Index py, Index pz, Index xi, Index eta, Index zeta, VisitorT&& visitor) {

				validateDiscretizationLimits(px, xi);
				validateDiscretizationLimits(py, eta);
				validateDiscretizationLimits(pz, zeta);

				using Traits = ElementTypeTraits<mesh::ElementFamily::Hex>;

				typename Traits::BasisType basisInst(px, py, pz);
				typename Traits::QuadratureVolumeType quadVolInst(xi, eta, zeta);

				const Index nBdy = (xi > eta) ? ((xi > zeta) ? xi : zeta) : ((eta > zeta) ? eta : zeta);
				typename Traits::QuadratureBoundaryType quadBdyInst(nBdy, nBdy);

				// propagates the visitor's return value, matching dispatchQuad() above
				return visitor.template operator()<NSD, 3, mesh::ElementFamily::Hex>(basisInst, quadVolInst, quadBdyInst);

			}

			template<Index NSD, typename VisitorT>
			bool dispatchNPD2(mesh::ElementFamily family, Index px, Index py, Index xi, Index eta, VisitorT&& visitor) {

				if (family == mesh::ElementFamily::Quad) {
					return dispatchQuad<NSD>(px, py, xi, eta, std::forward<VisitorT>(visitor));
				}

				return false;

			}

			template<Index NSD, typename VisitorT>
			bool dispatchNPD3(mesh::ElementFamily family, Index px, Index py, Index pz, Index xi, Index eta, Index zeta, VisitorT&& visitor) {

				if (family == mesh::ElementFamily::Hex) {
					return dispatchHex<NSD>(px, py, pz, xi, eta, zeta, std::forward<VisitorT>(visitor));
				}

				return false;

			}

			template<typename VisitorT>
			bool dispatch(Index nsd, Index npd, mesh::ElementFamily family, Index px, Index py, Index pz, Index xi, Index eta, Index zeta, VisitorT&& visitor) {

				if (npd == 2) {

					if (nsd == 2) {
						return dispatchNPD2<2>(family, px, py, xi, eta, std::forward<VisitorT>(visitor));
					}
					if (nsd == 3) {
						return dispatchNPD2<3>(family, px, py, xi, eta, std::forward<VisitorT>(visitor));
					}
					return false;

				}

				if (npd == 3) {

					if (nsd == 3) {
						return dispatchNPD3<3>(family, px, py, pz, xi, eta, zeta, std::forward<VisitorT>(visitor));
					}
					return false;

				}

				return false;

			}

		} // namespace dispatch
	} // namespace fem
} // namespace residuum

#endif
