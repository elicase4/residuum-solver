namespace residuum::fem::assembly {

	template<Index numDOFs, Index SpatialDim, GatherMode Mode, typename VectorT>
	PDE_HOST PDE_DEVICE void gatherElementVector(const Index* nodeIDs, Index nodesPerElement, const Real* nodeCoords, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry* bcRegistry, const Real time, Real* Ue, const VectorT* U){

		if constexpr (Mode == GatherMode::Free) {

			for (Index i = 0; i < nodesPerElement; ++i){
				for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){

					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					if (topoDOF.isConstrained(TdofIDi)) continue;
					Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);

					Ue[i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j] = U->data()[AdofIDi];

				}
			}

		} else {

			for (Index i = 0; i < nodesPerElement; ++i){

				Real bcVal[topology::TopologicalDOF<numDOFs>::dofsPerNode];
				bool haveBcVal = false;

				for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){

					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					
					// gather free nodes
					if (!topoDOF.isConstrained(TdofIDi)) {
						if constexpr (Mode == GatherMode::Full) {
							Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);
							Ue[i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j] = U->data()[AdofIDi];
						}
						continue;
					}

					// gather constrained nodes
					if (!haveBcVal) {
						Int rngTag = topoDOF.getConstraintTag(TdofIDi);
						const auto* entries = bcRegistry->getEntries(rngTag);
						if (entries) {
							for (const auto& entry : *entries) {
								entry->eval(time, &nodeCoords[SpatialDim*i], bcVal);
							}
						}
						haveBcVal = true;
					}

					Ue[i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j] = bcVal[j];

				}

			}

		}

	}

	template<Index numDOFs, ScatterMode Mode, typename VectorT>
	PDE_HOST PDE_DEVICE void scatterElementVector(const Index* nodeIDs, Index nodesPerElement, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real* Xe, VectorT& X, Real coefficient){

		for (Index i = 0; i < nodesPerElement; ++i){
			for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){

				Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
				if (topoDOF.isConstrained(TdofIDi)) continue;
				Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);

				X.data()[AdofIDi] += coefficient * Xe[i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j];

			}
		}

	}

	template<Index numDOFs, ScatterMode Mode, typename VectorT>
	PDE_HOST PDE_DEVICE void scatterElementDiagonal(const Index* nodeIDs, Index nodesPerElement, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real* Ke, VectorT& diag){

		constexpr Index dofsPerNode = topology::TopologicalDOF<numDOFs>::dofsPerNode;
		const Index localSize = nodesPerElement * dofsPerNode;

		for (Index i = 0; i < nodesPerElement; ++i){
			for (Index j = 0; j < dofsPerNode; ++j){

				Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
				if (topoDOF.isConstrained(TdofIDi)) continue;
				Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);
				Index localRow = i*dofsPerNode + j;

				diag.data()[AdofIDi] += Ke[localRow*localSize + localRow];

			}
		}

	}

	template<Index numDOFs, Index MaxNodesPerElement, ScatterMode Mode, typename MatrixT>
	PDE_HOST PDE_DEVICE void scatterElementMatrix(const Index* nodeIDs, Index nodesPerElement, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real* Ke, MatrixT& K, Real coefficient){

		constexpr Index dofsPerNode = topology::TopologicalDOF<numDOFs>::dofsPerNode;
		constexpr Index maxLocalSize = MaxNodesPerElement * numDOFs;
		const Index localSize = nodesPerElement * dofsPerNode;

		// this element's local->global column map
		Index globalCol[maxLocalSize];
		for (Index k = 0; k < nodesPerElement; ++k){
			for (Index l = 0; l < dofsPerNode; ++l){
				Index TdofIDk = topoDOF.getNodeDOF(nodeIDs[k], l);
				globalCol[k*dofsPerNode + l] = topoDOF.isConstrained(TdofIDk) ? Index(-1) : topoDOF.toAlgebraic(TdofIDk);
			}
		}

		// compact to the free entries
		Index sortedLocalCol[maxLocalSize];
		Index sortedGlobalCol[maxLocalSize];
		Index numFreeCols = 0;
		for (Index c = 0; c < localSize; ++c){
			if (globalCol[c] == Index(-1)) continue;
			sortedLocalCol[numFreeCols] = c;
			sortedGlobalCol[numFreeCols] = globalCol[c];
			++numFreeCols;
		}

		for (Index a = 1; a < numFreeCols; ++a){
			Index gcol = sortedGlobalCol[a];
			Index lcol = sortedLocalCol[a];
			Index b = a;
			while (b > 0 && sortedGlobalCol[b-1] > gcol){
				sortedGlobalCol[b] = sortedGlobalCol[b-1];
				sortedLocalCol[b] = sortedLocalCol[b-1];
				--b;
			}
			sortedGlobalCol[b] = gcol;
			sortedLocalCol[b] = lcol;
		}

		// merge this element's sorted target columns against K's already-sorted colIdx for that row in a single linear pass
		const Index* colIdx = K.colIdx();
		Real* data = K.data();
		const Index* rowPtr = K.rowPtr();

		for (Index i = 0; i < nodesPerElement; ++i){
			for (Index j = 0; j < dofsPerNode; ++j){

				Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
				if (topoDOF.isConstrained(TdofIDi)) continue;
				Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);
				Index localRow = i*dofsPerNode + j;

				Index rowEnd = rowPtr[AdofIDi + 1];
				Index p = rowPtr[AdofIDi];

				for (Index c = 0; c < numFreeCols; ++c){

					Index targetCol = sortedGlobalCol[c];
					while (p < rowEnd && colIdx[p] < targetCol) ++p;

					// colIdx is guaranteed to contain targetCol
					data[p] += coefficient * Ke[localRow*localSize + sortedLocalCol[c]];

				}

			}
		}

	}

} // namespace residuum::fem::assembly
