#ifndef RESIDUUM_SOLVER_CONFIG_EVALUATIONMETHOD_HPP
#define RESIDUUM_SOLVER_CONFIG_EVALUATIONMETHOD_HPP

namespace residuum {
	namespace solver {
		namespace config {

			// pointwise/strong (boundary quadrature) vs consistent/weak (volume residual reaction) evaluation of a boundary quantity; applies to any equation, not just heateq's flux
			enum class EvaluationMethod {
				Strong,
				Weak
			}; // enum class EvaluationMethod

		} // namespace config
	} // namespace solver
} // namespace residuum

#endif
