#ifndef RESIDUUM_EQUATION_HEATEQ_EVAL_TEMPERATUREFIELD_HPP
#define RESIDUUM_EQUATION_HEATEQ_EVAL_TEMPERATUREFIELD_HPP

#include <type_traits>

#include "fem/evaluator/EvalField.hpp"

namespace residuum::equation::heateq::evaluator {

	struct TemperatureField {

		static constexpr Index NumComponents = 1;

		PDE_HOST PDE_DEVICE void eval(const auto& qp, const Real* Ue, Real* outValue) const {

			Real u = 0;

			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				u += qp.N[a] * Ue[a];
			}

			outValue[0] = u;
		}

		PDE_HOST PDE_DEVICE void evalGradient(const auto& qp, const Real* Ue, Real* outGrad) const {

			using QP = std::decay_t<decltype(qp)>;

			for (Index i = 0; i < QP::SpatialDim; ++i){
				outGrad[i] = 0;
			}

			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				for (Index i = 0; i < QP::SpatialDim; ++i){
					outGrad[i] += qp.dNdx[a*QP::SpatialDim + i] * Ue[a];
				}
			}

		}

	}; // struct TemperatureField

} // namespace residuum::equation::heateq::evaluator

#endif
