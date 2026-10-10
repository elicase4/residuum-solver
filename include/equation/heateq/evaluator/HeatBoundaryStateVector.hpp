#ifndef RESIDUUM_EQUATION_HEATEQ_EVAL_HEATBOUNDARYSTATEVECTOR_HPP
#define RESIDUUM_EQUATION_HEATEQ_EVAL_HEATBOUNDARYSTATEVECTOR_HPP

#include <cstring>

#include "equation/heateq/evaluator/TemperatureField.hpp"
#include "fem/assembly/ElementMap.hpp"
#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/evaluator/FieldState.hpp"

namespace residuum {
	namespace equation {
		namespace heateq {
			namespace evaluator {

				template<Index numDOFs, Index SpatialDim, Index ParametricDim, typename VectorT>
				struct HeatBoundaryStateVector {

					const VectorT* valueSource = nullptr;

					Real Ue_value[fem::dispatch::kMaxNodesPerElement<ParametricDim> * numDOFs];

					// accessor
					PDE_HOST PDE_DEVICE Real* valueBuffer() { return Ue_value; }

					template<typename TopoDOFT>
					PDE_HOST PDE_DEVICE void gatherElement(const Index* nodeIDs, Index nodesPerElement, const Real* nodeCoords, const TopoDOFT& topoDOF, const fem::boundary::EssentialBoundaryRegistry* bcRegistry, Real time) {

						std::memset(Ue_value, 0, sizeof(Ue_value));
						fem::assembly::gatherElementVector<numDOFs, SpatialDim, fem::assembly::GatherMode::Full>(nodeIDs, nodesPerElement, nodeCoords, topoDOF, bcRegistry, time, Ue_value, valueSource);

					}

					template<typename QuadraturePointT>
					PDE_HOST PDE_DEVICE void interpolate(QuadraturePointT& qp) const {
						heateq::evaluator::TemperatureField().eval(qp, Ue_value, &qp.T.valueData[0]);
					}

				}; // struct HeatBoundaryStateVector

			} // namespace evaluator
		} // namespace heateq
	} // namespace equation
} // namespace residuum

#endif
