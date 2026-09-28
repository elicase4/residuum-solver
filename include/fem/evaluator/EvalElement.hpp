#ifndef RESIDUUM_FEM_EVAL_EVALELEMENT_HPP
#define RESIDUUM_FEM_EVAL_EVALELEMENT_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace residuum {
	namespace fem {
		namespace evaluator {

			template<typename ElementT>
			concept EvalElement = requires(ElementT elem, const Real* nodeCoords, const Real time) {

				{ ElementT::SpatialDim } -> std::convertible_to<Index>;
				{ ElementT::ParametricDim } -> std::convertible_to<Index>;

				{ elem.nodesPerElement() } -> std::convertible_to<Index>;
				{ elem.bindElement(nodeCoords, time) } -> std::same_as<void>;

			}; // concept EvalElement

		} // namespace evaluator
	} // namespace fem
} // namespace residuum

#endif
