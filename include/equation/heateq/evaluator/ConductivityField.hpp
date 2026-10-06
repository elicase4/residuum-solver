#ifndef RESIDUUM_EQUATION_HEATEQ_EVAL_CONDUCTIVITYFIELD_HPP
#define RESIDUUM_EQUATION_HEATEQ_EVAL_CONDUCTIVITYFIELD_HPP

#include "fem/evaluator/EvalField.hpp"

namespace residuum::equation::heateq::evaluator {

	// one representative scalar (not the full tensor, for anisotropic cases); assumes qp.K is
	// already populated (by ConductivityModel::eval(qp), same ordering forms already assume)
	struct ConductivityField {

		static constexpr Index NumComponents = 1;

		PDE_HOST PDE_DEVICE void eval(const auto& qp, const Real*, Real* outValue) const {
			outValue[0] = qp.K[0];
		}

		PDE_HOST PDE_DEVICE void evalGradient(const auto&, const Real*, Real*) const {} // unused

	}; // struct ConductivityField

} // namespace residuum::equation::heateq::evaluator

#endif
