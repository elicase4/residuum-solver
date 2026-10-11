#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include "fem/assembly/Assembler.hpp"

namespace residuum::fem::quantity {

template<>
class QuantityEvaluator<linalg::types::backend::CPU> {
public:

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointBoundary EvalQPT, typename ModelT, typename QuantityFormsT, typename QuadratureT, typename StateVectorT>
	requires evaluator::EvalModel<ModelT, EvalQPT>
	static void evaluateBoundaryRegistry(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const QuantityFormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, StateVectorT& state, BoundaryQuantityRegistry<QuantityFormsT::TotalComponents>& registry){

		registry.reset();

		EvalEleT localEle = evalEle;

		// boundary quadrature
		Real xi[fem::dispatch::kMaxQuadraturePointsTotalBoundary<EvalEleT::ParametricDim>*(EvalEleT::ParametricDim-1)];
		Real w[fem::dispatch::kMaxQuadraturePointsTotalBoundary<EvalEleT::ParametricDim>];
		assembly::materializeQuadrature(quadrature, xi, w);

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e) {

			// get element rngTags, skip elements with no face on ANY registered tag at all
			const Int* rngTags = mesh.getBoundaryTag(e);
			if (!registry.touchesAny(rngTags, mesh.data.facesPerElement)) continue;

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEleT::SpatialDim*fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];
			assembly::gatherElementCoords<EvalEleT::SpatialDim>(mesh, nodeIDs, localEle.nodesPerElement(), nodeCoords);

			state.gatherElement(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, &bcRegistry, time);

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// each face accumulates into its own tag's slot, if registered
			for (Index f = 0; f < mesh.data.facesPerElement; ++f) {

				if (!registry.hasTag(rngTags[f])) continue;

				auto& entry = registry.entries().at(rngTags[f]);

				// qp data
				EvalQPT qp(localEle, f);

				// quadrature loop
				for (Index q = 0; q < quadrature.numPointsTotal(); ++q){

					qp.evaluate(&xi[(EvalEleT::ParametricDim-1)*q], w[q]);
					state.interpolate(qp);
					model.eval(qp);
					model.evalGradient(qp);

					forms.computeElementLevelValue(qp, state.valueBuffer(), entry.result.data());

					accumulateMeasure<EvalEleT::SpatialDim>(qp, entry.measureSum);

				}

			}

		}

		// finalize Average-mode slices per tag, using each tag's OWN measureSum
		for (auto& [tag, entry] : registry.entries()) {
			forms.finalize(entry.result.data(), entry.measureSum);
		}

	}

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPVolT, evaluator::EvalQuadraturePointBoundary EvalQPBdyT, typename ModelT, typename QuadratureVolumeT, typename QuadratureBoundaryT, typename StateVectorT, typename... FormsWithOperandT>
	requires evaluator::EvalModel<ModelT, EvalQPVolT>
	static void evaluateBoundaryReaction(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const EvalEleT& evalEle, const QuadratureVolumeT& quadratureVolume, const QuadratureBoundaryT& quadratureBoundary, StateVectorT& state, BoundaryQuantityRegistry<numDOFs>& registry, FormsWithOperandT... pairs){

		registry.reset();

		// allocate Oe on the stack -- shared accumulator for every pair
		Real Oe[fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];

		EvalEleT localEle = evalEle;

		// which tag (if any) has already claimed each local node this element
		Int nodeClaimedByTag[fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];

		// volume quadrature
		Real xiVol[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>*(EvalEleT::ParametricDim)];
		Real wVol[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>];
		assembly::materializeQuadrature(quadratureVolume, xiVol, wVol);

		// boundary quadrature
		Real xiBdy[fem::dispatch::kMaxQuadraturePointsTotalBoundary<EvalEleT::ParametricDim>*(EvalEleT::ParametricDim-1)];
		Real wBdy[fem::dispatch::kMaxQuadraturePointsTotalBoundary<EvalEleT::ParametricDim>];
		assembly::materializeQuadrature(quadratureBoundary, xiBdy, wBdy);

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e) {

			// get element rngTags, skip elements with no face on ANY registered tag at all
			const Int* rngTags = mesh.getBoundaryTag(e);
			if (!registry.touchesAny(rngTags, mesh.data.facesPerElement)) continue;

			// zero-out Oe
			std::memset(Oe, 0.0, sizeof(Oe));

			// -1 means no tag has claimed this local node yet this element
			std::fill(nodeClaimedByTag, nodeClaimedByTag + localEle.nodesPerElement(), Int(-1));

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEleT::SpatialDim*fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];
			assembly::gatherElementCoords<EvalEleT::SpatialDim>(mesh, nodeIDs, localEle.nodesPerElement(), nodeCoords);

			state.gatherElement(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, &bcRegistry, time);

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// volume qp data
			EvalQPVolT qpVol(localEle);

			// volume quadrature loop
			for (Index q = 0; q < quadratureVolume.numPointsTotal(); ++q){

				qpVol.evaluate(&xiVol[EvalEleT::ParametricDim*q], wVol[q]);
				state.interpolate(qpVol);
				model.eval(qpVol);
				model.evalGradient(qpVol);

				// each pair reads its own selected buffer off state, accumulating into the same Oe
				(pairs.forms.computeElementLevelVector(qpVol, pairs.selector(state), Oe), ...);

			}

			// face loop
			for (Index f = 0; f < mesh.data.facesPerElement; ++f) {

				if (!registry.hasTag(rngTags[f])) continue;

				auto& entry = registry.entries().at(rngTags[f]);

				// qp data
				EvalQPBdyT qpBdy(localEle, f);

				for (Index q = 0; q < quadratureBoundary.numPointsTotal(); ++q){

					qpBdy.evaluate(&xiBdy[(EvalEleT::ParametricDim-1)*q], wBdy[q]);

					accumulateMeasure<EvalEleT::SpatialDim>(qpBdy, entry.measureSum);

				}

				for (Index i = 0; i < qpBdy.nodesPerFace(); ++i) {
					if (nodeClaimedByTag[qpBdy.faceNodeLocalIDs[i]] != rngTags[f]){
						for (Index c = 0; c < numDOFs; ++c) {
							entry.result.data()[c] += Oe[qpBdy.faceNodeLocalIDs[i]*numDOFs + c];
						}
						nodeClaimedByTag[qpBdy.faceNodeLocalIDs[i]] = rngTags[f];
					}
				}

			}

		}

	}

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FieldT, typename QuadratureT, typename RecoveryT, typename StateVectorT>
	requires evaluator::EvalModel<ModelT, EvalQPT> && recovery::NodalRecoveryStrategy<RecoveryT>
	static void evaluateNodalField(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const FieldT& field, const EvalEleT& evalEle, const QuadratureT& quadrature, StateVectorT& state, const RecoveryT& recovery, Real* nodalOut){

		constexpr Index NumComponents = FieldT::NumComponents;

		std::vector<Real> nodalWeight(mesh.data.numNodes, 0.0);
		std::fill(nodalOut, nodalOut + mesh.data.numNodes * NumComponents, 0.0);

		EvalEleT localEle = evalEle;

		Real xi[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>*EvalEleT::ParametricDim];
		Real w[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>];
		assembly::materializeQuadrature(quadrature, xi, w);

		// element loop -- every element, not just boundary-touching ones
		for (Index e = 0; e < mesh.data.numElements; ++e) {

			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEleT::SpatialDim*fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];
			assembly::gatherElementCoords<EvalEleT::SpatialDim>(mesh, nodeIDs, localEle.nodesPerElement(), nodeCoords);

			state.gatherElement(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, &bcRegistry, time);

			localEle.bindElement(nodeCoords, time);

			EvalQPT qp(localEle);

			for (Index q = 0; q < quadrature.numPointsTotal(); ++q){

				qp.evaluate(&xi[EvalEleT::ParametricDim*q], w[q]);
				state.interpolate(qp);
				model.eval(qp);
				model.evalGradient(qp);

				Real contribution[NumComponents];
				field.eval(qp, state.valueBuffer(), contribution);

				// scatter this quadrature point's contribution to each of its element's local nodes, weighted by that node's own shape function value here
				for (Index a = 0; a < localEle.nodesPerElement(); ++a) {
					const Index globalNode = nodeIDs[a];
					const Real weight = qp.N[a] * qp.measure * qp.w;
					recovery.accumulate(contribution, weight, &nodalOut[globalNode*NumComponents], &nodalWeight[globalNode], NumComponents);
				}

			}

		}

		for (Index n = 0; n < mesh.data.numNodes; ++n) {
			recovery.finalize(&nodalOut[n*NumComponents], nodalWeight[n], NumComponents);
		}

	}

}; // class QuantityEvaluator<linalg::types::backend::CPU>

} // namespace residuum::fem::quantity
