#ifndef RESIDUUM_EQUATION_HEATEQ_EVAL_EVALQUADRATUREPOINTVOLUME_HPP
#define RESIDUUM_EQUATION_HEATEQ_EVAL_EVALQUADRATUREPOINTVOLUME_HPP

#include "equation/heateq/evaluator/EvalField.hpp"
#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/evaluator/EvalQuadraturePointVolume.hpp"

namespace residuum::equation::heateq::evaluator {

	template<typename ElementT, typename BasisT, typename GeometryT>
	class EvalQuadraturePointVolume {
	public:

		ElementT element;

		EvalQuadraturePointVolume(const ElementT& elem) : element(elem) {}

		static constexpr Index SpatialDim = ElementT::SpatialDim;
		static constexpr Index ParametricDim = ElementT::ParametricDim;

		Index nodesPerElement() const { return element.nodesPerElement(); }

		// parent element attributes
		const Real time = element.t;
		const Real* coords = element.nodeCoords;

		// physical coordinate
		Real x[SpatialDim];

		// quadrature
		Real xi[ParametricDim];
		Real w;

		Real N[fem::dispatch::kMaxNodesPerElement<ParametricDim>];

		// ref gradients
		Real dNdxi[ParametricDim*fem::dispatch::kMaxNodesPerElement<ParametricDim>];

		// physical gradients
		Real dNdx[SpatialDim*fem::dispatch::kMaxNodesPerElement<ParametricDim>];

		// geometry
		Real J[SpatialDim*ParametricDim];
		Real g[ParametricDim*ParametricDim];

		// measure
		Real measure;

		// conductivity model
		Real K[SpatialDim*SpatialDim];
		Real dKdT[SpatialDim*SpatialDim];

		// heat capacity model
		Real rho;
		Real cp;
		Real dcpdT;

		// field state definition
		struct DOFFieldState {
			Real value;
			Real gradient[SpatialDim];
			Real rate;
		};

		// temperature field state
		DOFFieldState T;

		// number of auxiliary DOF states 
		static constexpr Index NumAuxStates = 1;

		PDE_HOST PDE_DEVICE void evaluate(const Real* xi_q, const Real weight){

			// set quad info
			for (Index pD = 0; pD < ParametricDim; ++pD){
				xi[pD] = xi_q[pD];
			}
			w = weight;

			// basis function evaluation
			element.basis().eval(xi, N);
			element.basis().evalGradient(xi, dNdxi);

			// geometry
			GeometryT::mapToPhysical(coords, N, x, nodesPerElement());
			GeometryT::computeJacobian(coords, dNdxi, J, nodesPerElement());
			GeometryT::computeMetric(J, g);
			measure = GeometryT::computeMeasure(g);

			// transforms
			GeometryT::transformGradient(J, g, dNdxi, dNdx, nodesPerElement());

		}

		PDE_HOST PDE_DEVICE void interpolateFields(const Real* Ue, const Real* const* auxStates) {
			EvalField().eval(*this, Ue, &T.value);
			EvalField().evalGradient(*this, Ue, T.gradient);
			EvalField().eval(*this, auxStates[0], &T.rate);
		}

	}; // class EvalQuadraturePointVolume

} // namespace residuum::equation::heateq::evaluator

#endif
