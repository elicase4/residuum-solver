#ifndef RESIDUUM_EQUATION_HEATEQ_FORM_MASSFORM_HPP
#define RESIDUUM_EQUATION_HEATEQ_FORM_MASSFORM_HPP

#include "fem/form/BilinearForm.hpp"

namespace residuum::equation::heateq {

	template<typename QuadraturePointVolumeT>
	struct MassForm {

		PDE_HOST PDE_DEVICE static void computeElementLevelMatrix(const QuadraturePointVolumeT& qp, Real* Me){
			
			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				for (Index b = 0; b < qp.nodesPerElement(); ++b){
					Me[a * qp.nodesPerElement() + b] += qp.rho * qp.cp * qp.N[a] * qp.N[b] * qp.measure * qp.w;
				}
			}
		}

		PDE_HOST PDE_DEVICE static void computeElementLevelVector(const QuadraturePointVolumeT& qp, const Real* Ue, Real* Oe){

			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				for (Index b = 0; b < qp.nodesPerElement(); ++b){
					Oe[a] += qp.rho * qp.cp * qp.N[a] * qp.N[b] * Ue[b] * qp.measure * qp.w;
				}
			}
		}

	}; // struct MassForm

} // namespace residuum::equation::heateq

#endif
