namespace residuum::topology {

template<Index numDOFs>
TopologicalDOF<numDOFs>::TopologicalDOF(const mesh::Mesh& mesh, fem::dof::DOFOrdering ordering) : mesh_(mesh), ordering_(ordering), numFreeDOFsPerField_(dofsPerNode) {
	
	numGlobalDOFs_ = mesh_.data.numNodes * dofsPerNode;
	numFreeDOFs_ = numGlobalDOFs_;
	
	// initialize mapping
	topoToAlg_.resize(numGlobalDOFs_);
	for (Index i = 0; i < numGlobalDOFs_; ++i){
		topoToAlg_[i] = static_cast<Int>(i);
	}

}

template<Index numDOFs>
void TopologicalDOF<numDOFs>::getElementDOFs(Index elemId, Index* dofs) const {

	const Index* nodes = mesh_.getElementNodes(elemId);
	const Index npe = mesh_.data.nodesPerElement;

	Index k = 0;
	for (Index a = 0; a < npe; ++a) {
		for (Index c = 0; c < dofsPerNode; ++c) {
			dofs[k++] = nodes[a] * dofsPerNode + c;
		}
	}

}

template<Index numDOFs>
template<typename BasisT>
void TopologicalDOF<numDOFs>::buildConstraints(const BasisT& basis, const fem::boundary::EssentialBoundaryRegistry& bcRegistry){

	std::set<Index> constrainedSet;
	Index faceNodes[fem::dispatch::kMaxNodesPerElement<BasisT::ParametricDim>];

	// loop over elements and faces to mark constrained dofs
	for (Index e = 0; e < mesh_.data.numElements; ++e) {
		for (Index f = 0; f < mesh_.data.facesPerElement; ++f){
			
			// check if element has boundary nodes
			if (!mesh_.isOnBoundary(e,f)) continue;

			// check for essential bcs
			const Int* tagPtr = mesh_.getBoundaryTag(e);
			Int tag = tagPtr[f];

			// get face nodes
			Index npf = basis.nodesPerFace(f);
			basis.getFaceNodes(f, faceNodes);

			// get element nodes
			const Index* elemNodes = mesh_.getElementNodes(e);

			// mark all constrained DOFs on face f
			for (Index i = 0; i < npf; ++i) {
				
				Index globalNode = elemNodes[faceNodes[i]];
				
				for (Index c = 0; c < dofsPerNode; ++c){
					
					if (!bcRegistry.isEssential(tag, c)) continue;

					Index dof = getNodeDOF(globalNode, c);
					
					constrainedSet.insert(dof);
					constraintTags_[dof] = tag;
				}
			
			}

		}
	}

	// build algebraic numering from free dofs
	algToTopo_.clear();
	algToTopo_.reserve(numGlobalDOFs_ - constrainedSet.size());
	
	for (Index i = 0; i < numGlobalDOFs_; ++i) topoToAlg_[i] = -1;
	Index algIndex = 0;
	
	if (ordering_ == fem::dof::DOFOrdering::Interleaved) {
		
		for (Index topoDOF = 0; topoDOF < numGlobalDOFs_; ++topoDOF) {
			if (constrainedSet.count(topoDOF)) continue;
			topoToAlg_[topoDOF] = static_cast<Int>(algIndex); // map to free dof
			algToTopo_.push_back(topoDOF);
			++algIndex;
		}

	} else {

		for (Index c = 0; c < dofsPerNode; ++c) {
			for (Index node = 0; node < mesh_.data.numNodes; ++node) {
				Index topoDOF = getNodeDOF(node, c);
				if (constrainedSet.count(topoDOF)) continue;
				topoToAlg_[topoDOF] = static_cast<Int>(algIndex); // map to free dof
				algToTopo_.push_back(topoDOF);
				++algIndex;
			}
		}

	}

	// set num of free dofs to final tally
	numFreeDOFs_ = algIndex;
	numFreeDOFsPerField_ = (dofsPerNode > 0) ? (numFreeDOFs_ / dofsPerNode ) : 0; // extension later to mixed order elements

}

template<Index numDOFs>
Int TopologicalDOF<numDOFs>::getConstraintTag(Index topoDOF) const {
	auto it = constraintTags_.find(topoDOF);
	return (it != constraintTags_.end()) ? it->second : -1;
}

} // namespace residuum::topology
