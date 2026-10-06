#ifndef RESIDUUM_EQUATION_HEATEQ_EVAL_HEATFLUXFIELD_HPP
#define RESIDUUM_EQUATION_HEATEQ_EVAL_HEATFLUXFIELD_HPP

#include "equation/heateq/evaluator/TemperatureField.hpp"
#include "fem/evaluator/EvalField.hpp"

namespace residuum::equation::heateq::evaluator {

	// Q = -K*gradT; assumes qp.K is already populated (by ConductivityModel::eval(qp), same
	// ordering every other qp.K reader in this codebase already assumes)
	template<Index SpatialDim>
	struct HeatFluxField {

		static constexpr Index NumComponents = SpatialDim;

		PDE_HOST PDE_DEVICE void eval(const auto& qp, const Real* Ue, Real* outValue) const {

			Real gradT[SpatialDim];
			TemperatureField().evalGradient(qp, Ue, gradT);

			for (Index i = 0; i < SpatialDim; ++i) {
				outValue[i] = 0;
				for (Index j = 0; j < SpatialDim; ++j) {
					outValue[i] -= qp.K[i*SpatialDim + j] * gradT[j];
				}
			}

		}

		PDE_HOST PDE_DEVICE void evalGradient(const auto&, const Real*, Real*) const {} // unused

	}; // struct HeatFluxField

} // namespace residuum::equation::heateq::evaluator

#endif
