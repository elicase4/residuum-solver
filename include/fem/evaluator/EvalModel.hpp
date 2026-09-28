#ifndef RESIDUUM_FEM_EVAL_EVALMODEL_HPP
#define RESIDUUM_FEM_EVAL_EVALMODEL_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace residuum {
	namespace fem {
		namespace evaluator {

			template<typename ModelT, typename QuadraturePointT>
			concept EvalModel = requires (const ModelT m, QuadraturePointT& qp) {

				{ m.eval(qp) } -> std::same_as<void>;
				{ m.evalGradient(qp) } -> std::same_as<void>;

			}; // concept EvalModel

		} // namespace evaluator
	} // namespace fem
} // namespace residuum

#endif
