#ifndef RESIDUUM_FEM_EVALUATOR_STATEVECTORSELECTORS_HPP
#define RESIDUUM_FEM_EVALUATOR_STATEVECTORSELECTORS_HPP

#include "config/Platform.hpp"
#include "core/Types.hpp"

namespace residuum {
	namespace fem {
		namespace evaluator {

			// FormsWithOperand selectors: named, device-callable in place of ad-hoc lambdas, which
			// lambdas aren't by default without their own PDE_HOST/PDE_DEVICE annotation. Generic
			// over any StateVectorT exposing the named buffer, so equation-agnostic code (solver::stage,
			// solver::problem) can use these without depending on a specific equation module.
			struct SelectValue {
				template<typename StateVectorT>
				PDE_HOST PDE_DEVICE Real* operator()(StateVectorT& s) const { return s.valueBuffer(); }
			}; // struct SelectValue

			struct SelectRate {
				template<typename StateVectorT>
				PDE_HOST PDE_DEVICE Real* operator()(StateVectorT& s) const { return s.rateBuffer(); }
			}; // struct SelectRate

			struct SelectOperand {
				template<typename StateVectorT>
				PDE_HOST PDE_DEVICE Real* operator()(StateVectorT& s) const { return s.operandBuffer(); }
			}; // struct SelectOperand

		} // namespace evaluator
	} // namespace fem
} // namespace residuum

#endif
