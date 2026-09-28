#ifndef RESIDUUM_EQUATION_HEATEQ_FORM_DIFFUSIONFORM_HPP
#define RESIDUUM_EQUATION_HEATEQ_FORM_DIFFUSIONFORM_HPP

#include "fem/form/BilinearForm.hpp"

namespace residuum::equation::heateq {

	template<typename QuadraturePointVolumeT>
	struct DiffusionForm {

		PDE_HOST PDE_DEVICE static void computeElementLevelMatrix(const QuadraturePointVolumeT& qp, Real* Ke){
			
			Real integrand, matvecprod;
			
			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				for (Index b = 0; b < qp.nodesPerElement(); ++b){
					
					integrand = 0.0;
					for (Index sDi = 0; sDi < qp.SpatialDim; ++sDi){
						matvecprod = 0.0;
						for (Index sDj = 0; sDj < qp.SpatialDim; ++sDj){
							matvecprod += qp.K[sDi*qp.SpatialDim + sDj] * qp.dNdx[b*qp.SpatialDim + sDj];
						}
						integrand += qp.dNdx[a*qp.SpatialDim + sDi] * matvecprod;
					}

					Ke[a * qp.nodesPerElement() + b] += integrand * qp.measure * qp.w;
				}
			}
		
		}
		
		PDE_HOST PDE_DEVICE static void computeElementLevelVector(const QuadraturePointVolumeT& qp, const Real* Ue, Real* Oe){
			
			Real integrand, matvecprod;

			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				for (Index b = 0; b < qp.nodesPerElement(); ++b){
					
					integrand = 0.0;
					for (Index sDi = 0; sDi < qp.SpatialDim; ++sDi){
						matvecprod = 0.0;
						for (Index sDj = 0; sDj < qp.SpatialDim; ++sDj){
							matvecprod += qp.K[sDi*qp.SpatialDim + sDj] * qp.dNdx[b*qp.SpatialDim + sDj];
						}
						integrand += qp.dNdx[a*qp.SpatialDim + sDi] * matvecprod;
					}
					
					Oe[a] += integrand * Ue[b] * qp.measure * qp.w;
				}
			}
		}

	}; // struct DiffusionForm

} // namespace residuum::equation::heateq

#endif
