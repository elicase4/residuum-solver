#ifndef RESIDUUM_FEM_EVAL_EVALFUNCTION_HPP
#define RESIDUUM_FEM_EVAL_EVALFUNCTION_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include <concepts>

namespace residuum {
	namespace fem {
		namespace evaluator {

			template<typename FunctionT>
			concept EvalFunction = requires (const FunctionT f, const Real time, const Real* x, Real* outValue, Real* outGrad) {
				
				{ FunctionT::NumComponents } -> std::convertible_to<Index>;
				{ FunctionT::SpatialDim } -> std::convertible_to<Index>;
				{ f.eval(time, x, outValue) } -> std::same_as<void>;
				{ f.evalGradient(time, x, outGrad) } -> std::same_as<void>;
			
			}; // concept EvalFunction

		} // namespace evaluator
	} // namespace fem
} // namespace residuum

#endif
