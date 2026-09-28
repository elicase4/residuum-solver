#ifndef RESIDUUM_FEM_EVAL_EVALNODALDATA_HPP
#define RESIDUUM_FEM_EVAL_EVALNODALDATA_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace residuum {
	namespace fem {
		namespace evaluator {

			// node-indexed, no interpolation; distinct from EvalField, which interpolates a quadrature-point value
			template<typename NodalDataT>
			concept EvalNodalData = requires (const NodalDataT d, Index nodeID, Real* outValue) {

				{ NodalDataT::NumComponents } -> std::convertible_to<Index>;
				{ d.eval(nodeID, outValue) } -> std::same_as<void>;

			}; // concept EvalNodalData

		} // namespace evaluator
	} // namespace fem
} // namespace residuum

#endif
