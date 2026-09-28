namespace residuum::fem::boundary {

template<>
class BoundaryApplicator<linalg::types::backend::CPU> {
public:

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FormsT, typename QuadratureT>
	requires evaluator::EvalModel<ModelT, EvalQPT>
	static void applyEssentialBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const EssentialBoundaryRegistry& bcRegistry, const Real time, const ModelT& model, const FormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

		// allocate Fe on the stack
		Real Fe[(fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*topology::TopologicalDOF<numDOFs>::dofsPerNode)];

		// allocate Ge on the stack
		Real Ge[(fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*topology::TopologicalDOF<numDOFs>::dofsPerNode)];

		EvalEleT localEle = evalEle;

		// quadrature points/weights
		Real xi[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>*EvalEleT::ParametricDim];
		Real w[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>];
		quadrature.getPoints(xi);
		quadrature.getWeights(w);

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){

			// zero-out Fe
			std::memset(Fe, 0.0, sizeof(Fe));

			// zero-out Ge
			std::memset(Ge, 0.0, sizeof(Ge));

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEleT::SpatialDim * fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEleT::SpatialDim; ++sD){
					nodeCoords[EvalEleT::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// qp data
			EvalQPT qp(localEle);

			// fill Ge
			fem::assembly::gatherElementVector<numDOFs, EvalEleT::SpatialDim, fem::assembly::GatherMode::Constrained>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, &bcRegistry, time, Ge);

			// gather any form-specific element data
			forms.gatherElementData(nodeIDs, localEle.nodesPerElement());

			// quadrature loop
			for (Index q = 0; q < quadrature.numPointsTotal(); ++q){
				qp.evaluate(&xi[EvalEleT::ParametricDim*q], w[q]);
				model.eval(qp);
				model.evalGradient(qp);
				forms.computeElementLevelVector(qp, Ge, Fe);
			}

			// scatter Fe into F
			fem::assembly::scatterElementVector<numDOFs, fem::assembly::ScatterMode::Free>(nodeIDs, localEle.nodesPerElement(), topoDOF, Fe, F, Real(-1));

		}

	}

	template<Index numDOFs, evaluator::EvalElement EvalEleT, typename EvalQPT, typename QuadratureT>
	static void applyNaturalBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const NaturalBoundaryRegistry<EvalQPT>& bcRegistry, const Real time, const EvalEleT& evalEle, const QuadratureT& quadrature, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

		// allocate local space for Fe
		Real Fe[(fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*topology::TopologicalDOF<numDOFs>::dofsPerNode)];

		EvalEleT localEle = evalEle;

		// boundary quadrature
		Real xi[fem::dispatch::kMaxQuadraturePointsTotalBoundary<EvalEleT::ParametricDim>*(EvalEleT::ParametricDim-1)];
		Real w[fem::dispatch::kMaxQuadraturePointsTotalBoundary<EvalEleT::ParametricDim>];
		quadrature.getPoints(xi);
		quadrature.getWeights(w);

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEleT::SpatialDim*fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEleT::SpatialDim; ++sD){
					nodeCoords[EvalEleT::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// get element rngTags
			const Int* rngTags = mesh.getBoundaryTag(e);

			// face loop
			for (Index f = 0; f < mesh.data.facesPerElement; ++f){

				// zero-out Fe
				std::memset(Fe, 0.0, sizeof(Fe));

				// get face rng tag
				Int rngTag = rngTags[f];
				if (rngTag < 0) continue;
				if (!bcRegistry.hasAny(rngTag)) continue;

				// get face nodes information
				const Index nodesPerFace = localEle.basis().nodesPerFace(f);
				
				Index faceNodeLocalIDs[fem::dispatch::kMaxNodesPerElementBoundary<EvalEleT::ParametricDim>];
				Index faceNodeGlobalIDs[fem::dispatch::kMaxNodesPerElementBoundary<EvalEleT::ParametricDim>];

				localEle.basis().getFaceNodes(f, faceNodeLocalIDs);
				const Index* elemNodeGlobalIDs = mesh.getElementNodes(e);

				for (Index i = 0; i < nodesPerFace; ++i){
					faceNodeGlobalIDs[i] = elemNodeGlobalIDs[faceNodeLocalIDs[i]];
				}

				// gather any operator-specific face data
				bcRegistry.gatherFaceElementData(rngTag, faceNodeGlobalIDs, nodesPerFace);

				// qp data
				EvalQPT qp(localEle, f);

				// quadrature loop
				for (Index q = 0; q < quadrature.numPointsTotal(); ++q){

					qp.evaluate(&xi[(EvalEleT::ParametricDim-1)*q], w[q]);
					bcRegistry.apply(rngTag, qp, Fe);

				}

				// scatter Fe into F
				fem::assembly::scatterElementVector<numDOFs, fem::assembly::ScatterMode::Free>(faceNodeGlobalIDs, nodesPerFace, topoDOF, Fe, F);

			}

		}


	}

}; // class BoundaryApplicator

} // namespace residuum::fem::boundary
