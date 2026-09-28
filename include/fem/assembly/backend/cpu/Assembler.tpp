namespace residuum::fem::assembly {

template<>
class Assembler<linalg::types::backend::CPU> {
public:

	template<Index numDOFs>
	static linalg::types::CSRMatrix<Real, linalg::types::backend::CPU> createMatrix(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF){

		constexpr Index dofsPerNode = topology::TopologicalDOF<numDOFs>::dofsPerNode;
		const Index nRows = topoDOF.numFreeDOFs();

		linalg::types::CSRMatrix<Real, linalg::types::backend::CPU> K(nRows, nRows);

		// visits every (row, col) sparsity touch an element contributes
		auto forEachTouch = [&](auto&& visit){

			for (Index e = 0; e < mesh.data.numElements; ++e){

				const Index* nodeIDs = mesh.getElementNodes(e);

				for (Index i = 0; i < mesh.data.nodesPerElement; ++i){
					for (Index j = 0; j < dofsPerNode; ++j){

						Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
						if (topoDOF.isConstrained(TdofIDi)) continue;
						Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);

						for (Index k = 0; k < mesh.data.nodesPerElement; ++k){
							for (Index l = 0; l < dofsPerNode; ++l){

								Index TdofIDk = topoDOF.getNodeDOF(nodeIDs[k], l);
								if (topoDOF.isConstrained(TdofIDk)) continue;
								Index AdofIDk = topoDOF.toAlgebraic(TdofIDk);

								visit(AdofIDi, AdofIDk);

							}
						}

					}
				}

			}

		};

		// pass 1: count raw (pre-dedup) touches per row
		std::vector<Index> rowOffset(nRows + 1, 0);
		forEachTouch([&](Index row, Index){ ++rowOffset[row + 1]; });
		for (Index i = 0; i < nRows; ++i) rowOffset[i + 1] += rowOffset[i];

		// pass 2: fill one flat buffer via a per-row cursor
		std::vector<Index> rawAdj(rowOffset[nRows]);
		std::vector<Index> cursor(rowOffset.begin(), rowOffset.end() - 1);
		forEachTouch([&](Index row, Index col){ rawAdj[cursor[row]++] = col; });

		// dedupe each row in place (sort + unique on a slice of the one flat buffer
		K.rowPtr()[0] = 0;
		for (Index i = 0; i < nRows; ++i){
			Index* begin = rawAdj.data() + rowOffset[i];
			Index* end = rawAdj.data() + rowOffset[i + 1];
			std::sort(begin, end);
			Index* uniqueEnd = std::unique(begin, end);
			K.rowPtr()[i + 1] = K.rowPtr()[i] + static_cast<Index>(uniqueEnd - begin);
		}

		// allocate K's real (exactly-sized) storage, then copy each row's deduplicated entries
		K.resize(K.rowPtr()[nRows]);
		for (Index i = 0; i < nRows; ++i){
			Index count = K.rowPtr()[i + 1] - K.rowPtr()[i];
			std::copy(rawAdj.data() + rowOffset[i], rawAdj.data() + rowOffset[i] + count, K.colIdx() + K.rowPtr()[i]);
		}

		return K;
	}

	template<Index numDOFs>
	static linalg::types::Vector<Real, linalg::types::backend::CPU> createVector(const mesh::Mesh&, const topology::TopologicalDOF<numDOFs>& topoDOF){

		linalg::types::Vector<Real, linalg::types::backend::CPU> F(topoDOF.numFreeDOFs());
		return F;

	}

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FormsT, typename QuadratureT, GatherMode Mode>
	requires evaluator::EvalModel<ModelT, EvalQPT>
	static void assembleMatrix(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const ModelT& model, const FormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, linalg::types::backend::CPU>& U, const std::array<const linalg::types::Vector<Real, linalg::types::backend::CPU>*, EvalQPT::NumAuxStates>& auxStates, linalg::types::CSRMatrix<Real, linalg::types::backend::CPU>& K, const fem::boundary::EssentialBoundaryRegistry* bcRegistry){

		// allocate Ke on the stack
		Real Ke[(fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs) * (fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs)];

		// allocate Ue on the stack
		Real Ue[(fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs)];

		// allocate the per-element auxiliary-state buffers on the stack
		Real Ue_aux[EvalQPT::NumAuxStates][(fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs)];
		const Real* Ue_auxPtrs[EvalQPT::NumAuxStates];

		// zero-out data in K
		K.zero();

		EvalEleT localEle = evalEle;

		// setup quadrature points and weights
		Real xi[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>*EvalEleT::ParametricDim];
		Real w[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>];
		quadrature.getPoints(xi);
		quadrature.getWeights(w);

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){

			// zero-out Ke
			std::memset(Ke, 0.0, sizeof(Ke));

			// zero-out Ue
			std::memset(Ue, 0.0, sizeof(Ue));

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEleT::SpatialDim * fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEleT::SpatialDim; ++sD){
					nodeCoords[EvalEleT::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// gather U into Ue
			gatherElementVector<numDOFs, EvalEleT::SpatialDim, Mode>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, bcRegistry, time, Ue, &U);

			// aux states are rate fields, so a constrained node contributes 0 rather than a looked-up value
			for (Index s = 0; s < EvalQPT::NumAuxStates; ++s) {
				std::memset(Ue_aux[s], 0.0, sizeof(Ue_aux[s]));
				if (auxStates[s] != nullptr) {
					gatherElementVector<numDOFs, EvalEleT::SpatialDim, GatherMode::Free>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, bcRegistry, time, Ue_aux[s], auxStates[s]);
				}
				Ue_auxPtrs[s] = Ue_aux[s];
			}

			// gather any form-specific element data
			forms.gatherElementData(nodeIDs, localEle.nodesPerElement());

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
				forms.computeElementLevelMatrix(qp, Ke);
			}

			// scatter Ke into K
			scatterElementMatrix<numDOFs, fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>, ScatterMode::Free>(nodeIDs, localEle.nodesPerElement(), topoDOF, Ke, K);

		}

	}

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FormsT, typename QuadratureT, GatherMode Mode>
	requires evaluator::EvalModel<ModelT, EvalQPT>
	static void assembleDiagonal(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const ModelT& model, const FormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, linalg::types::backend::CPU>* fieldSource, const std::array<const linalg::types::Vector<Real, linalg::types::backend::CPU>*, EvalQPT::NumAuxStates>& auxStates, linalg::types::Vector<Real, linalg::types::backend::CPU>& diag, const fem::boundary::EssentialBoundaryRegistry* bcRegistry){

		// allocate Ke on the stack -- the full local matrix is computed (reusing
		// computeElementLevelMatrix as-is, no new form method needed), only its diagonal is kept
		Real Ke[(fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs) * (fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs)];

		// allocate Ue on the stack
		Real Ue[(fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs)];

		// allocate the per-element auxiliary-state buffers on the stack
		Real Ue_aux[EvalQPT::NumAuxStates][(fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs)];
		const Real* Ue_auxPtrs[EvalQPT::NumAuxStates];

		// zero-out diag
		diag.zero();

		EvalEleT localEle = evalEle;

		// setup quadrature points and weights
		Real xi[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>*EvalEleT::ParametricDim];
		Real w[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>];
		quadrature.getPoints(xi);
		quadrature.getWeights(w);

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){

			// zero-out Ke
			std::memset(Ke, 0.0, sizeof(Ke));

			// zero-out Ue
			std::memset(Ue, 0.0, sizeof(Ue));

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEleT::SpatialDim * fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEleT::SpatialDim; ++sD){
					nodeCoords[EvalEleT::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// a null fieldSource means the model doesn't depend on the field at all
			// (e.g. a constant-property steady solve), so Ue is left zeroed
			if (fieldSource != nullptr) {
				gatherElementVector<numDOFs, EvalEleT::SpatialDim, Mode>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, bcRegistry, time, Ue, fieldSource);
			}

			// aux states are rate fields, so a constrained node contributes 0 rather than a looked-up value
			for (Index s = 0; s < EvalQPT::NumAuxStates; ++s) {
				std::memset(Ue_aux[s], 0.0, sizeof(Ue_aux[s]));
				if (auxStates[s] != nullptr) {
					gatherElementVector<numDOFs, EvalEleT::SpatialDim, GatherMode::Free>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, bcRegistry, time, Ue_aux[s], auxStates[s]);
				}
				Ue_auxPtrs[s] = Ue_aux[s];
			}

			// gather any form-specific element data
			forms.gatherElementData(nodeIDs, localEle.nodesPerElement());

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
				forms.computeElementLevelMatrix(qp, Ke);
			}

			// scatter only Ke's diagonal into diag -- no sparsity structure needed
			scatterElementDiagonal<numDOFs, ScatterMode::Free>(nodeIDs, localEle.nodesPerElement(), topoDOF, Ke, diag);

		}

	}

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FormsT, typename QuadratureT, GatherMode Mode>
	requires evaluator::EvalModel<ModelT, EvalQPT>
	static void assembleVector(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const ModelT& model, const FormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, linalg::types::backend::CPU>& U, const linalg::types::Vector<Real, linalg::types::backend::CPU>* fieldSource, const std::array<const linalg::types::Vector<Real, linalg::types::backend::CPU>*, EvalQPT::NumAuxStates>& auxStates, linalg::types::Vector<Real, linalg::types::backend::CPU>& F, const fem::boundary::EssentialBoundaryRegistry* bcRegistry){

		// allocate Fe on the stack
		Real Fe[fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];

		// allocate Ue on the stack
		Real Ue[fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];

		// allocate a buffer for fieldSource, used only when fieldSource != nullptr
		Real Ue_lin[fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];

		// allocate the per-element auxiliary-state buffers on the stack
		Real Ue_aux[EvalQPT::NumAuxStates][fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];
		const Real* Ue_auxPtrs[EvalQPT::NumAuxStates];

		// zero-out data in F
		F.zero();

		// local mutable copy
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

			// zero-out Ue
			std::memset(Ue, 0.0, sizeof(Ue));

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEleT::SpatialDim * fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEleT::SpatialDim; ++sD){
					nodeCoords[EvalEleT::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// gather U into Ue
			gatherElementVector<numDOFs, EvalEleT::SpatialDim, Mode>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, bcRegistry, time, Ue, &U);

			// gather fieldSource when it differs from the operand
			const Real* fieldPtr = Ue;
			if (fieldSource != nullptr) {
				std::memset(Ue_lin, 0.0, sizeof(Ue_lin));
				gatherElementVector<numDOFs, EvalEleT::SpatialDim, GatherMode::Full>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, bcRegistry, time, Ue_lin, fieldSource);
				fieldPtr = Ue_lin;
			}

			// aux states are rate fields, so a constrained node contributes 0 rather than a looked-up value
			for (Index s = 0; s < EvalQPT::NumAuxStates; ++s) {
				std::memset(Ue_aux[s], 0.0, sizeof(Ue_aux[s]));
				if (auxStates[s] != nullptr) {
					gatherElementVector<numDOFs, EvalEleT::SpatialDim, GatherMode::Free>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, bcRegistry, time, Ue_aux[s], auxStates[s]);
				}
				Ue_auxPtrs[s] = Ue_aux[s];
			}

			// gather any form-specific element data
			forms.gatherElementData(nodeIDs, localEle.nodesPerElement());

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// qp data
			EvalQPT qp(localEle);

			// quadrature loop
			for (Index q = 0; q < quadrature.numPointsTotal(); ++q){
				qp.evaluate(&xi[EvalEleT::ParametricDim*q], w[q]);
				qp.interpolateFields(fieldPtr, Ue_auxPtrs);
				model.eval(qp);
				model.evalGradient(qp);
				forms.computeElementLevelVector(qp, Ue, Fe);
			}

			// scatter Fe into F
			scatterElementVector<numDOFs, ScatterMode::Free>(nodeIDs, localEle.nodesPerElement(), topoDOF, Fe, F);

		}

	}

}; // class Assembler <linalg::types::backend::CPU>

} // namespace residuum::fem::assembly
