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

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FormsT, typename QuadratureT, typename StateVectorT>
	requires evaluator::EvalModel<ModelT, EvalQPT> && evaluator::StateVector<StateVectorT, numDOFs, EvalQPT>
	static void assembleMatrix(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const ModelT& model, const FormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, StateVectorT& state, linalg::types::CSRMatrix<Real, linalg::types::backend::CPU>& K, const fem::boundary::EssentialBoundaryRegistry* bcRegistry){

		// allocate Ke on the stack
		Real Ke[(fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs) * (fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs)];

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

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEleT::SpatialDim * fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];
			gatherElementCoords<EvalEleT::SpatialDim>(mesh, nodeIDs, localEle.nodesPerElement(), nodeCoords);

			state.gatherElement(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, bcRegistry, time);

			// gather any form-specific element data
			forms.gatherElementData(nodeIDs, localEle.nodesPerElement());

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// qp data
			EvalQPT qp(localEle);

			// quadrature loop
			for (Index q = 0; q < quadrature.numPointsTotal(); ++q){
				qp.evaluate(&xi[EvalEleT::ParametricDim*q], w[q]);
				state.interpolate(qp);
				model.eval(qp);
				model.evalGradient(qp);
				forms.computeElementLevelMatrix(qp, Ke);
			}

			// scatter Ke into K
			scatterElementMatrix<numDOFs, fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>, ScatterMode::Free>(nodeIDs, localEle.nodesPerElement(), topoDOF, Ke, K);

		}

	}

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FormsT, typename QuadratureT, typename StateVectorT>
	requires evaluator::EvalModel<ModelT, EvalQPT> && evaluator::StateVector<StateVectorT, numDOFs, EvalQPT>
	static void assembleDiagonal(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const ModelT& model, const FormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, StateVectorT& state, linalg::types::Vector<Real, linalg::types::backend::CPU>& diag, const fem::boundary::EssentialBoundaryRegistry* bcRegistry){

		// allocate Ke on the stack -- the full local matrix is computed (reusing
		// computeElementLevelMatrix as-is, no new form method needed), only its diagonal is kept
		Real Ke[(fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs) * (fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs)];

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

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEleT::SpatialDim * fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];
			gatherElementCoords<EvalEleT::SpatialDim>(mesh, nodeIDs, localEle.nodesPerElement(), nodeCoords);

			state.gatherElement(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, bcRegistry, time);

			// gather any form-specific element data
			forms.gatherElementData(nodeIDs, localEle.nodesPerElement());

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// qp data
			EvalQPT qp(localEle);

			// quadrature loop
			for (Index q = 0; q < quadrature.numPointsTotal(); ++q){
				qp.evaluate(&xi[EvalEleT::ParametricDim*q], w[q]);
				state.interpolate(qp);
				model.eval(qp);
				model.evalGradient(qp);
				forms.computeElementLevelMatrix(qp, Ke);
			}

			// scatter only Ke's diagonal into diag -- no sparsity structure needed
			scatterElementDiagonal<numDOFs, ScatterMode::Free>(nodeIDs, localEle.nodesPerElement(), topoDOF, Ke, diag);

		}

	}

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename QuadratureT, typename StateVectorT, typename... FormsWithOperandT>
	requires evaluator::EvalModel<ModelT, EvalQPT> && evaluator::StateVector<StateVectorT, numDOFs, EvalQPT>
	static void assembleVector(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const ModelT& model, const EvalEleT& evalEle, const QuadratureT& quadrature, StateVectorT& state, linalg::types::Vector<Real, linalg::types::backend::CPU>& F, const fem::boundary::EssentialBoundaryRegistry* bcRegistry, FormsWithOperandT... pairs){

		// allocate Fe on the stack -- shared accumulator for every pair
		Real Fe[fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];

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

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEleT::SpatialDim * fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];
			gatherElementCoords<EvalEleT::SpatialDim>(mesh, nodeIDs, localEle.nodesPerElement(), nodeCoords);

			// gather and interpolate state once, shared by every pair
			state.gatherElement(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, bcRegistry, time);

			// gather any form-specific element data, for every pair's forms
			(pairs.forms.gatherElementData(nodeIDs, localEle.nodesPerElement()), ...);

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// qp data
			EvalQPT qp(localEle);

			// quadrature loop
			for (Index q = 0; q < quadrature.numPointsTotal(); ++q){
				qp.evaluate(&xi[EvalEleT::ParametricDim*q], w[q]);
				state.interpolate(qp);
				model.eval(qp);
				model.evalGradient(qp);
				(pairs.forms.computeElementLevelVector(qp, pairs.selector(state), Fe), ...);
			}

			// scatter the combined Fe into F
			scatterElementVector<numDOFs, ScatterMode::Free>(nodeIDs, localEle.nodesPerElement(), topoDOF, Fe, F);

		}

	}

}; // class Assembler <linalg::types::backend::CPU>

} // namespace residuum::fem::assembly
