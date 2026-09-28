#ifndef RESIDUUM_FEM_EVAL_EVALFIELD_HPP
#define RESIDUUM_FEM_EVAL_EVALFIELD_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace residuum {
	namespace fem {
		namespace evaluator {

			template<typename FieldT, typename QuadraturePointT>
			concept EvalField = requires (const FieldT field, QuadraturePointT& qp, const Real* Ue, Real* outValue, Real* outGrad) {

				{ FieldT::NumComponents } -> std::convertible_to<Index>;

				{ field.eval(qp, Ue, outValue) } -> std::same_as<void>;
				{ field.evalGradient(qp, Ue, outGrad) } -> std::same_as<void>;

			}; // concept EvalField

		} // namespace evaluator
	} // namespace fem
} // namespace residuum

#endif
