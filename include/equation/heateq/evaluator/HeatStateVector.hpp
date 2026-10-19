#ifndef RESIDUUM_EQUATION_HEATEQ_EVAL_HEATSTATEVECTOR_HPP
#define RESIDUUM_EQUATION_HEATEQ_EVAL_HEATSTATEVECTOR_HPP

#include <cstring>

#include "equation/heateq/evaluator/TemperatureField.hpp"
#include "fem/assembly/ElementMap.hpp"
#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/evaluator/FieldState.hpp"
#include "fem/evaluator/StateVectorSelectors.hpp"

namespace residuum {
	namespace equation {
		namespace heateq {
			namespace evaluator {

				template<Index numDOFs, Index SpatialDim, Index ParametricDim, typename VectorT>
				struct HeatStateVector {

					const VectorT* valueSource = nullptr;
					const VectorT* rateSource = nullptr;
					const VectorT* operandSource = nullptr;

					Real Ue_value[fem::dispatch::kMaxNodesPerElement<ParametricDim> * numDOFs];
					Real Ue_rate[fem::dispatch::kMaxNodesPerElement<ParametricDim> * numDOFs];
					Real Ue_operand[fem::dispatch::kMaxNodesPerElement<ParametricDim> * numDOFs];

					// accessors
					PDE_HOST PDE_DEVICE Real* valueBuffer() { return Ue_value; }
					PDE_HOST PDE_DEVICE Real* rateBuffer() { return Ue_rate; }
					PDE_HOST PDE_DEVICE Real* operandBuffer() { return Ue_operand; }

					// constrainedValues is the flat, per-call-lowered Dirichlet buffer (see fem::boundary::lowerEssentialBoundaryValues) -- gatherElementVector's Full-mode path reads it directly, no registry lookup here
					template<typename TopoDOFT>
					PDE_HOST PDE_DEVICE void gatherElement(const Index* nodeIDs, Index nodesPerElement, const Real* nodeCoords, const TopoDOFT& topoDOF, const Real* constrainedValues, Real time) {

						(void) nodeCoords; (void) time;

						std::memset(Ue_value, 0, sizeof(Ue_value));
						fem::assembly::gatherElementVector<numDOFs, fem::assembly::GatherMode::Full>(nodeIDs, nodesPerElement, topoDOF, Ue_value, valueSource, constrainedValues);

						std::memset(Ue_rate, 0, sizeof(Ue_rate));
						if (rateSource != nullptr) {
							fem::assembly::gatherElementVector<numDOFs, fem::assembly::GatherMode::Free>(nodeIDs, nodesPerElement, topoDOF, Ue_rate, rateSource);
						}

						// a Krylov basis vector is free-dof-only by construction -- gathering it never needs Full/Constrained handling
						std::memset(Ue_operand, 0, sizeof(Ue_operand));
						if (operandSource != nullptr) {
							fem::assembly::gatherElementVector<numDOFs, fem::assembly::GatherMode::Free>(nodeIDs, nodesPerElement, topoDOF, Ue_operand, operandSource);
						}

					}

					template<typename QuadraturePointT>
					PDE_HOST PDE_DEVICE void interpolate(QuadraturePointT& qp) const {
						heateq::evaluator::TemperatureField().eval(qp, Ue_value, &qp.T.valueData[0]);
						heateq::evaluator::TemperatureField().evalGradient(qp, Ue_value, qp.gradT);
						heateq::evaluator::TemperatureField().eval(qp, Ue_rate, &qp.T.derivatives[0][0]);
					}

				}; // struct HeatStateVector

			} // namespace evaluator
		} // namespace heateq
	} // namespace equation
} // namespace residuum

#endif
