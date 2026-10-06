#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include "fem/assembly/Assembler.hpp"

namespace residuum::fem::quantity {

template<>
class QuantityEvaluator<linalg::types::backend::CPU> {
public:

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename QuantityFormsT, typename QuadratureT>
	requires evaluator::EvalModel<ModelT, EvalQPT>
	static void evaluateDomain(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const QuantityFormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, linalg::types::backend::CPU>& U, const std::array<const linalg::types::Vector<Real, linalg::types::backend::CPU>*, EvalQPT::NumAuxStates>& auxStates, Real* out){

		// zero-out output buffer & measure summation
		for (Index c = 0; c < QuantityFormsT::TotalComponents; ++c) out[c] = 0.0;
		Real measureSum = 0.0;

		// allocate Ue on the stack
		Real Ue[fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];

		// allocate the per-element auxiliary-state buffers on the stack
		Real Ue_aux[EvalQPT::NumAuxStates][fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];
		const Real* Ue_auxPtrs[EvalQPT::NumAuxStates];

		EvalEleT localEle = evalEle;

		// quadrature points/weights
		Real xi[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>*EvalEleT::ParametricDim];
		Real w[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>];
		quadrature.getPoints(xi);
		quadrature.getWeights(w);

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e) {

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEleT::SpatialDim * fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEleT::SpatialDim; ++sD){
					nodeCoords[EvalEleT::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// gather the complete nodal solution
			fem::assembly::gatherElementVector<numDOFs, EvalEleT::SpatialDim, fem::assembly::GatherMode::Full>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, &bcRegistry, time, Ue, &U);

			// aux states are rate fields, so a constrained node contributes 0 rather than a looked-up value
			for (Index s = 0; s < EvalQPT::NumAuxStates; ++s) {
				std::memset(Ue_aux[s], 0.0, sizeof(Ue_aux[s]));
				if (auxStates[s] != nullptr) {
					fem::assembly::gatherElementVector<numDOFs, EvalEleT::SpatialDim, fem::assembly::GatherMode::Free>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, &bcRegistry, time, Ue_aux[s], auxStates[s]);
				}
				Ue_auxPtrs[s] = Ue_aux[s];
			}

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// qp data
			EvalQPT qp(localEle);

			// quadrature loop
			for (Index q = 0; q < quadrature.numPointsTotal(); ++q){

				qp.evaluate(&xi[EvalEleT::ParametricDim*q], w[q]);
				qp.interpolateFields(Ue, Ue_auxPtrs);
				model.eval(qp);
				model.evalGradient(qp);

				forms.computeElementLevelValue(qp, Ue, out);

				measureSum += qp.measure * qp.w;

			}

		}

		// finalize Average-mode slices
		forms.finalize(out, measureSum);

	}

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointBoundary EvalQPT, typename ModelT, typename QuantityFormsT, typename QuadratureT>
	requires evaluator::EvalModel<ModelT, EvalQPT>
	static void evaluateBoundaryRegistry(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const QuantityFormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, linalg::types::backend::CPU>& U, BoundaryQuantityRegistry<QuantityFormsT::TotalComponents>& registry){

		// zero every registered tag's slot & measureSum
		for (auto& [tag, entry] : registry.entries()) {
			entry.result.fill(0.0);
			entry.measureSum = 0.0;
		}

		// allocate Ue on the stack
		Real Ue[fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];

		EvalEleT localEle = evalEle;

		// boundary quadrature
		Real xi[fem::dispatch::kMaxQuadraturePointsTotalBoundary<EvalEleT::ParametricDim>*(EvalEleT::ParametricDim-1)];
		Real w[fem::dispatch::kMaxQuadraturePointsTotalBoundary<EvalEleT::ParametricDim>];
		quadrature.getPoints(xi);
		quadrature.getWeights(w);

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e) {

			// get element rngTags, skip elements with no face on ANY registered tag at all
			const Int* rngTags = mesh.getBoundaryTag(e);

			bool touchesAny = false;
			for (Index f = 0; f < mesh.data.facesPerElement; ++f) {
				if (registry.hasTag(rngTags[f])) { touchesAny = true; break; }
			}
			if (!touchesAny) continue;

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEleT::SpatialDim*fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEleT::SpatialDim; ++sD){
					nodeCoords[EvalEleT::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			fem::assembly::gatherElementVector<numDOFs, EvalEleT::SpatialDim, fem::assembly::GatherMode::Full>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, &bcRegistry, time, Ue, &U);

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
					qp.interpolateFields(Ue);
					model.eval(qp);
					model.evalGradient(qp);

					forms.computeElementLevelValue(qp, Ue, entry.result.data());

					Real normalMag = 0;
					for (Index sD = 0; sD < EvalEleT::SpatialDim; ++sD) {
						normalMag += qp.normal[sD] * qp.normal[sD];
					}
					entry.measureSum += std::sqrt(normalMag) * qp.w;

				}

			}

		}

		// finalize Average-mode slices per tag, using each tag's OWN measureSum
		for (auto& [tag, entry] : registry.entries()) {
			forms.finalize(entry.result.data(), entry.measureSum);
		}

	}

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPVolT, evaluator::EvalQuadraturePointBoundary EvalQPBdyT, typename ModelT, typename FormsT, typename QuadratureVolumeT, typename  QuadratureBoundaryT>
	requires evaluator::EvalModel<ModelT, EvalQPVolT>
	static void evaluateBoundaryReaction(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const FormsT& forms, const EvalEleT& evalEle, const QuadratureVolumeT& quadratureVolume, const QuadratureBoundaryT& quadratureBoundary, const linalg::types::Vector<Real, linalg::types::backend::CPU>& gatherSource, const std::array<const linalg::types::Vector<Real, linalg::types::backend::CPU>*, EvalQPVolT::NumAuxStates>& auxStates, bool zeroRegistry, bool computeMeasure, BoundaryQuantityRegistry<numDOFs>& registry){

		// a later call (e.g. the mass pass, after the stiffness+source pass already ran) accumulates
		// on top of the first call's result instead of re-zeroing it
		if (zeroRegistry) {
			for (auto& [tag, entry] : registry.entries()) {
				entry.result.fill(0.0);
				entry.measureSum = 0.0;
			}
		}

		// allocate Ue on the stack
		Real Ue[fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];

		// allocate Oe on the stack
		Real Oe[fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];

		// allocate the per-element auxiliary-state buffers on the stack
		Real Ue_aux[EvalQPVolT::NumAuxStates][fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];
		const Real* Ue_auxPtrs[EvalQPVolT::NumAuxStates];

		EvalEleT localEle = evalEle;

		// which tag (if any) has already claimed each local node this element -- dedupes a node
		// shared by two faces of the SAME tag on this element, without suppressing the same node's
		// (independent) contribution to a DIFFERENT tag's registry entry
		Int nodeClaimedByTag[fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];

		// volume quadrature
		Real xiVol[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>*(EvalEleT::ParametricDim)];
		Real wVol[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>];
		quadratureVolume.getPoints(xiVol);
		quadratureVolume.getWeights(wVol);

		// boundary quadrature
		Real xiBdy[fem::dispatch::kMaxQuadraturePointsTotalBoundary<EvalEleT::ParametricDim>*(EvalEleT::ParametricDim-1)];
		Real wBdy[fem::dispatch::kMaxQuadraturePointsTotalBoundary<EvalEleT::ParametricDim>];
		quadratureBoundary.getPoints(xiBdy);
		quadratureBoundary.getWeights(wBdy);
		
		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e) {

			// get element rngTags, skip elements with no face on ANY registered tag at all
			const Int* rngTags = mesh.getBoundaryTag(e);

			bool touchesAny = false;
			for (Index f = 0; f < mesh.data.facesPerElement; ++f) {
				if (registry.hasTag(rngTags[f])) { touchesAny = true; break; }
			}
			if (!touchesAny) continue;

			// zero-out Oe
			std::memset(Oe, 0.0, sizeof(Oe));

			// -1 means no tag has claimed this local node yet this element
			std::fill(nodeClaimedByTag, nodeClaimedByTag + localEle.nodesPerElement(), Int(-1));

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEleT::SpatialDim*fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEleT::SpatialDim; ++sD){
					nodeCoords[EvalEleT::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			fem::assembly::gatherElementVector<numDOFs, EvalEleT::SpatialDim, fem::assembly::GatherMode::Full>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, &bcRegistry, time, Ue, &gatherSource);

			// aux states are rate fields, so a constrained node contributes 0 rather than a looked-up value
			for (Index s = 0; s < EvalQPVolT::NumAuxStates; ++s) {
				std::memset(Ue_aux[s], 0.0, sizeof(Ue_aux[s]));
				if (auxStates[s] != nullptr) {
					fem::assembly::gatherElementVector<numDOFs, EvalEleT::SpatialDim, fem::assembly::GatherMode::Free>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, &bcRegistry, time, Ue_aux[s], auxStates[s]);
				}
				Ue_auxPtrs[s] = Ue_aux[s];
			}

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// volume qp data
			EvalQPVolT qpVol(localEle);

			// volume quadrature loop
			for (Index q = 0; q < quadratureVolume.numPointsTotal(); ++q){

				qpVol.evaluate(&xiVol[EvalEleT::ParametricDim*q], wVol[q]);
				qpVol.interpolateFields(Ue, Ue_auxPtrs);
				model.eval(qpVol);
				model.evalGradient(qpVol);

				forms.computeElementLevelVector(qpVol, Ue, Oe);

			}

			// face loop
			for (Index f = 0; f < mesh.data.facesPerElement; ++f) {

				if (!registry.hasTag(rngTags[f])) continue;

				auto& entry = registry.entries().at(rngTags[f]);
				
				// qp data
				EvalQPBdyT qpBdy(localEle, f);

				// measureSum is pure boundary geometry, independent of which gather source/forms are
				// active this call, so only the first (computeMeasure) call needs to accumulate it
				if (computeMeasure) {
					for (Index q = 0; q < quadratureBoundary.numPointsTotal(); ++q){

						qpBdy.evaluate(&xiBdy[(EvalEleT::ParametricDim-1)*q], wBdy[q]);

						Real normalMag = 0;
						for (Index sD = 0; sD < EvalEleT::SpatialDim; ++sD) {
							normalMag += qpBdy.normal[sD] * qpBdy.normal[sD];
						}
						entry.measureSum += std::sqrt(normalMag) * qpBdy.w;

					}
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

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FieldT, typename QuadratureT, typename RecoveryT>
	requires evaluator::EvalModel<ModelT, EvalQPT> && recovery::NodalRecoveryStrategy<RecoveryT>
	static void evaluateNodalField(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const FieldT& field, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, linalg::types::backend::CPU>& U, const RecoveryT& recovery, Real* nodalOut){

		constexpr Index NumComponents = FieldT::NumComponents;

		std::vector<Real> nodalWeight(mesh.data.numNodes, 0.0);
		std::fill(nodalOut, nodalOut + mesh.data.numNodes * NumComponents, 0.0);

		// allocate Ue on the stack
		Real Ue[fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];

		// output fields don't need a real rate, so the aux-state buffer stays zeroed
		Real Ue_aux[EvalQPT::NumAuxStates][fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];
		const Real* Ue_auxPtrs[EvalQPT::NumAuxStates];
		for (Index s = 0; s < EvalQPT::NumAuxStates; ++s) {
			std::memset(Ue_aux[s], 0.0, sizeof(Ue_aux[s]));
			Ue_auxPtrs[s] = Ue_aux[s];
		}

		EvalEleT localEle = evalEle;

		Real xi[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>*EvalEleT::ParametricDim];
		Real w[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>];
		quadrature.getPoints(xi);
		quadrature.getWeights(w);

		// element loop -- every element, not just boundary-touching ones
		for (Index e = 0; e < mesh.data.numElements; ++e) {

			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEleT::SpatialDim*fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEleT::SpatialDim; ++sD){
					nodeCoords[EvalEleT::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			fem::assembly::gatherElementVector<numDOFs, EvalEleT::SpatialDim, fem::assembly::GatherMode::Full>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, &bcRegistry, time, Ue, &U);

			localEle.bindElement(nodeCoords, time);

			EvalQPT qp(localEle);

			for (Index q = 0; q < quadrature.numPointsTotal(); ++q){

				qp.evaluate(&xi[EvalEleT::ParametricDim*q], w[q]);
				qp.interpolateFields(Ue, Ue_auxPtrs);
				model.eval(qp);
				model.evalGradient(qp);

				Real contribution[NumComponents];
				field.eval(qp, Ue, contribution);

				// scatter this quadrature point's contribution to each of its element's local nodes,
				// weighted by that node's own shape function value here (a lumped projection)
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
