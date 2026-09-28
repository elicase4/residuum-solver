namespace residuum::fem::dof {

PDE_HOST PDE_DEVICE Index AlgebraicDOF::getElementDOFs(const Index* elemTopoDOFs, const Index numTopoDOFs, const Int* topoToAlg, Index* elemAlgDOFs){
	
	Index count = 0;
	
	for (Index i = 0; i < numTopoDOFs; ++i) {
		Int algDOF = topoToAlg[elemTopoDOFs[i]];
		if (algDOF != -1){
			elemAlgDOFs[count++] = (Index) algDOF;
		}
	}

	return count;

}

} // namespace residuum::fem::dof
