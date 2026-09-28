#ifndef RESIDUUM_FEM_EVAL_EVALQUADRATUREPOINTBOUNDARY_HPP
#define RESIDUUM_FEM_EVAL_EVALQUADRATUREPOINTBOUNDARY_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace residuum {
	namespace fem {
		namespace evaluator {

			template<typename QuadraturePointBoundaryT>
			concept EvalQuadraturePointBoundary = requires(QuadraturePointBoundaryT qp, const Real* xi, const Real w) {

				{ QuadraturePointBoundaryT::SpatialDim } -> std::convertible_to<Index>;
				{ QuadraturePointBoundaryT::ParametricDim } -> std::convertible_to<Index>;

				{ qp.nodesPerElement() } -> std::convertible_to<Index>;
				{ qp.nodesPerFace() } -> std::convertible_to<Index>;
				{ qp.evaluate(xi, w) } -> std::same_as<void>;

				// face-local node data, indexed a=0..nodesPerFace()-1
				{ qp.faceNodeLocalIDs[0] } -> std::convertible_to<Index>;
				{ qp.Nface[0] } -> std::convertible_to<Real>;

				{ qp.normal[0] } -> std::convertible_to<Real>;
				{ qp.x[0] } -> std::convertible_to<Real>;
				{ qp.time } -> std::convertible_to<Real>;
				{ qp.w } -> std::convertible_to<Real>;

			}; // concept EvalQuadraturePointBoundary

		} // namespace evaluator
	} // namespace fem
} // namespace residuum

#endif
