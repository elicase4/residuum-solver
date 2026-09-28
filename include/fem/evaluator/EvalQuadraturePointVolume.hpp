#ifndef RESIDUUM_FEM_EVAL_EVALQUADRATUREPOINTVOLUME_HPP
#define RESIDUUM_FEM_EVAL_EVALQUADRATUREPOINTVOLUME_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace residuum {
	namespace fem {
		namespace evaluator {

			template<typename QuadraturePointVolumeT>
			concept EvalQuadraturePointVolume = requires(QuadraturePointVolumeT qp, const Real* xi, const Real w) {

				{ QuadraturePointVolumeT::SpatialDim } -> std::convertible_to<Index>;
				{ QuadraturePointVolumeT::ParametricDim } -> std::convertible_to<Index>;

				{ qp.nodesPerElement() } -> std::convertible_to<Index>;
				{ qp.evaluate(xi, w) } -> std::same_as<void>;

			}; // concept EvalQuadraturePointVolume

		} // namespace evaluator
	} // namespace fem
} // namespace residuum

#endif
