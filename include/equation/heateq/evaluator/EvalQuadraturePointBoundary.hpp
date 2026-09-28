#ifndef RESIDUUM_EQUATION_HEATEQ_EVAL_EVALQUADRATUREPOINTBOUNDARY_HPP
#define RESIDUUM_EQUATION_HEATEQ_EVAL_EVALQUADRATUREPOINTBOUNDARY_HPP

#include "equation/heateq/evaluator/EvalField.hpp"
#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/evaluator/EvalQuadraturePointBoundary.hpp"

namespace residuum::equation::heateq::evaluator {

	template<typename ElementT, typename BasisT, typename GeometryT>
	class EvalQuadraturePointBoundary {
	public:

		static constexpr Index SpatialDim = ElementT::SpatialDim;
		static constexpr Index ParametricDim = ElementT::ParametricDim;

		ElementT element;
		Int faceID;

		Index faceNodeLocalIDs[fem::dispatch::kMaxNodesPerElementBoundary<ParametricDim>];

		EvalQuadraturePointBoundary(const ElementT& elem, const Int fID) : element(elem), faceID(fID) {

			element.basis().getFaceNodes(faceID, faceNodeLocalIDs);

		}

		Index nodesPerElement() const { return element.nodesPerElement(); }

		Index nodesPerFace() const { return element.basis().nodesPerFace(faceID); }

		// parent element attributes
		const Real time = element.t;
		const Real* coords = element.nodeCoords;

		// physical coordinates
		Real x[SpatialDim];

		// reference coordinate
		Real xi[ParametricDim];

		// quadrature
		Real xi_face[ParametricDim-1];
		Real w;

		Real N[fem::dispatch::kMaxNodesPerElement<ParametricDim>];

		Real Nface[fem::dispatch::kMaxNodesPerElementBoundary<ParametricDim>];

		// ref gradients
		Real dNdxi[ParametricDim*fem::dispatch::kMaxNodesPerElement<ParametricDim>];

		// normal vectors
		Real normal[SpatialDim];
		Real normalRef[ParametricDim];

		// geometry
		Real J[SpatialDim*ParametricDim];
		Real g[ParametricDim*ParametricDim];

		// physical gradient
		Real dNdx[SpatialDim*fem::dispatch::kMaxNodesPerElement<ParametricDim>];

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
		};

		// temperature field state
		DOFFieldState T;

		PDE_HOST PDE_DEVICE void evaluate(const Real* xi_face_q, const Real weight){

			// set quad info
			for (Index pD = 0; pD < (ParametricDim - 1); ++pD){
				xi_face[pD] = xi_face_q[pD];
			}
			w = weight;

			// get volume parametric coordinates
			element.basis().mapFaceToElement(faceID, xi_face, xi);
			element.basis().eval(xi, N);
			element.basis().evalGradient(xi, dNdxi);
			element.basis().getFaceTopology(faceID, normalRef);

			// gather N onto face-local indices
			for (Index a = 0; a < nodesPerFace(); ++a){
				Nface[a] = N[faceNodeLocalIDs[a]];
			}

			// geometry
			GeometryT::mapToPhysical(coords, N, x, nodesPerElement());
			GeometryT::computeJacobian(coords, dNdxi, J, nodesPerElement());
			GeometryT::computeBoundaryNormal(J, normalRef, normal);
			GeometryT::computeMetric(J, g);

			// transforms
			GeometryT::transformGradient(J, g, dNdxi, dNdx, nodesPerElement());

		}

		PDE_HOST PDE_DEVICE void interpolateFields(const Real* Ue) {
			EvalField().eval(*this, Ue, &T.value);
			EvalField().evalGradient(*this, Ue, T.gradient);
		}

	}; // class EvalQuadraturePointBoundary

} // namespace residuum::equation::heateq::evaluator

#endif
