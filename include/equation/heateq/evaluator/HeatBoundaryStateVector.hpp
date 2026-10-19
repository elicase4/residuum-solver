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

					// constrainedValues is the flat, per-call-lowered Dirichlet buffer (see fem::boundary::lowerEssentialBoundaryValues) -- gatherElementVector's Full-mode path reads it directly, no registry lookup here
					template<typename TopoDOFT>
					PDE_HOST PDE_DEVICE void gatherElement(const Index* nodeIDs, Index nodesPerElement, const Real* nodeCoords, const TopoDOFT& topoDOF, const Real* constrainedValues, Real time) {

						(void) nodeCoords; (void) time;

						std::memset(Ue_value, 0, sizeof(Ue_value));
						fem::assembly::gatherElementVector<numDOFs, fem::assembly::GatherMode::Full>(nodeIDs, nodesPerElement, topoDOF, Ue_value, valueSource, constrainedValues);

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
