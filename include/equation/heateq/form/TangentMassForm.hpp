#ifndef RESIDUUM_EQUATION_HEATEQ_FORM_TANGENTMASSFORM_HPP
#define RESIDUUM_EQUATION_HEATEQ_FORM_TANGENTMASSFORM_HPP

#include "fem/form/BilinearForm.hpp"

namespace residuum::equation::heateq {

	template<typename QuadraturePointVolumeT>
	struct TangentMassForm {

		PDE_HOST PDE_DEVICE static void computeElementLevelMatrix(const QuadraturePointVolumeT& qp, Real* Me){

			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				for (Index b = 0; b < qp.nodesPerElement(); ++b){
					Me[a * qp.nodesPerElement() + b] += qp.rho * qp.dcpdT * qp.T.rate * qp.N[a] * qp.N[b] * qp.measure * qp.w;
				}
			}
		}

		PDE_HOST PDE_DEVICE static void computeElementLevelVector(const QuadraturePointVolumeT& qp, const Real* Ue, Real* Oe){

			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				for (Index b = 0; b < qp.nodesPerElement(); ++b){
					Oe[a] += qp.rho * qp.dcpdT * qp.T.rate * qp.N[a] * qp.N[b] * Ue[b] * qp.measure * qp.w;
				}
			}
		}

	}; // struct TangentMassForm

} // namespace residuum::equation::heateq

#endif
