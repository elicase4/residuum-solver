#include <cmath>
#include <cstring>

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
	static void evaluateBoundaryRegistry(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const QuantityFormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, linalg::types::backend::CPU>& U, BoundaryQuantityRegistry<QuantityFormsT>& registry){

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

}; // class QuantityEvaluator<linalg::types::backend::CPU>

} // namespace residuum::fem::quantity
