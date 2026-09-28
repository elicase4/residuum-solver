#include "mesh/exchange/gmsh/MeshConverter.hpp"

void residuum::mesh::exchange::gmsh::MeshConverter::toSolverMesh(residuum::mesh::Mesh& mesh, const residuum::mesh::exchange::gmsh::IntermediateMesh& input, const std::unordered_map<Int, Int>& physicalGroupMap) {

	if (input.empty()){
		throw std::runtime_error("MeshConverter: input IntermediateMesh is empty");
	}

	mesh.clear();
	std::unordered_map<Index, Index> gmshToSolver;
	residuum::mesh::exchange::gmsh::ElementType cellType = buildConnectivity(mesh, input, gmshToSolver);
	buildBoundaryTags(mesh, input, cellType, gmshToSolver, physicalGroupMap);

	if (!mesh.isValid()){
		throw std::runtime_error("MeshConverter: produced invalid mesh");
	}

}

residuum::mesh::exchange::gmsh::ElementType residuum::mesh::exchange::gmsh::MeshConverter::buildConnectivity(residuum::mesh::Mesh& mesh, const residuum::mesh::exchange::gmsh::IntermediateMesh& input, std::unordered_map<Index, Index>& gmshToSolver) {

	const ElementBlock* protoBlock = nullptr;
	
	Index totalElems = 0;
	
	// compute the total number of solver elements
	for (const auto* eb : input.cellBlocks()) {
	
		// check that cell blocks are assigned and are of uniform type
		if (!protoBlock) {
			protoBlock = eb;
		} else if (eb->type != protoBlock->type) {
			throw std::runtime_error("MeshConverter: mixed cell types are not supported");
		}

		totalElems += eb->elementIDs.size();

	}

	// set mesh metadata
	mesh.data.parametricDim = input.parametricDim;
	mesh.data.spatialDim = input.spatialDim;
	mesh.data.basisOrder = residuum::mesh::exchange::gmsh::basisOrder(protoBlock->type);
	mesh.data.elementFamily = residuum::mesh::exchange::gmsh::elementFamily(protoBlock->type);
	mesh.data.basisType = residuum::mesh::exchange::gmsh::basisType(protoBlock->type);
	mesh.data.nodesPerElement = protoBlock->nodesPerElement;
	mesh.data.facesPerElement = residuum::mesh::exchange::gmsh::facesPerElement(protoBlock->type);

	// set coordinates
	const Index numNodes = input.xyz.size() / 3;
	
	mesh.data.numNodes = numNodes;
	mesh.data.xyz.resize(numNodes * input.spatialDim);

	gmshToSolver.reserve(numNodes);
	
	Index solverNode = 0;
	Index xyzOffset = 0;

	for (const auto& nb: input.nodeBlocks){

		const Index nBlockNodes = nb.nodeIDs.size();

		for (Index i = 0; i < nBlockNodes; ++i){

			const Index gmshTag = nb.nodeIDs[i];
			gmshToSolver[gmshTag] = solverNode;

			for (Index d = 0; d < input.spatialDim; ++d){
				mesh.data.xyz[solverNode * input.spatialDim + d] = input.xyz[3*xyzOffset + d];
			}

			++solverNode;
			++xyzOffset;

		}

	}

	// set connectivity
	mesh.data.numElements = totalElems;
	mesh.data.ien.resize(totalElems * protoBlock->nodesPerElement);
	Index elemOffset = 0;
	for (const auto* eb : input.cellBlocks()){

		const Index nElem = eb->elementIDs.size();
		for (Index e = 0; e < nElem; ++e){
			
			// get element connectivity
			const Index* conn = &eb->connectivity[e * eb->nodesPerElement];
			
			// reorder connectivity to row major
			auto reordered = reorderConnectivity(conn, eb->type);
			
			// store node id for each connectivity entry
			for (Index n = 0; n < eb->nodesPerElement; ++n){

				auto it = gmshToSolver.find(reordered[n]);
				if (it == gmshToSolver.end()) {
					throw std::runtime_error("MeshConverter: missing gmsh node tag");
				}

				mesh.data.ien[(elemOffset + e)*mesh.data.nodesPerElement + n] = it->second;

			}

		}

		elemOffset += nElem;

	}

	return protoBlock->type;
	
}

void residuum::mesh::exchange::gmsh::MeshConverter::buildBoundaryTags(residuum::mesh::Mesh& mesh, const residuum::mesh::exchange::gmsh::IntermediateMesh& input, residuum::mesh::exchange::gmsh::ElementType cellType, std::unordered_map<Index, Index>& gmshToSolver, const std::unordered_map<Int, Int>& physicalGroupMap) {

	const Index fpe = mesh.data.facesPerElement;
	const Index nElem = mesh.data.numElements;
	
	// Initialize rng to -1 before filling
	mesh.data.rng.assign(nElem * fpe, -1);

	// Build face-set -> (elemID, localFace) lookup from volume mesh. Key is sorted node IDs of a face
	std::unordered_map<std::vector<Index>, std::pair<Index, Index>, residuum::mesh::exchange::gmsh::VecHash> faceMap;
	faceMap.reserve(nElem * fpe * 2);

	// fill in the face map
	for (Index e = 0; e < nElem; ++e){
		
		const Index* elemNodes = mesh.getElementNodes(e);
		
		for (Index f = 0; f < fpe; ++f){
			auto faceNodes = residuum::mesh::exchange::gmsh::localFaceNodes(elemNodes, cellType, f);
			std::sort(faceNodes.begin(), faceNodes.end());
			faceMap[faceNodes] = {e, f};

		}
	
	}

	// iterate over the boundary element blocks
	for (const auto* eb : input.boundaryBlocks()) {

		// resolve solver boundary tag with input map
		Int solverTag = eb->physicalTag; // use raw physical tag by default
		
		// check input physicalGroupMap for solver boundary ID
		auto mapIt = physicalGroupMap.find(eb->physicalTag);
		if (mapIt != physicalGroupMap.end()){
			solverTag = mapIt->second;
		}
		
		// loop over face elements in this block
		const Index nFaceElem = eb->elementIDs.size();
		for (Index fe = 0; fe < nFaceElem; ++fe){

			// get connectivity key for each node and reorder
			const Index* conn = &eb->connectivity[fe*eb->nodesPerElement];
			std::vector<Index> key = reorderConnectivity(conn, eb->type);

			// map gmsh node-based key to solver node-based key
			for (Index i  = 0; i < key.size(); ++i){

				auto it = gmshToSolver.find(key[i]);
				
				if (it == gmshToSolver.end()) {
					throw std::runtime_error("MeshConverter: boundary node tag missing");
				}

				key[i] = it->second;

			}
			
			// sort output
			std::sort(key.begin(), key.end());
			
			// get info about face node from face map via gmsh face node id lookup
			auto it = faceMap.find(key);

			if (it == faceMap.end()){
				continue;
			}

			auto [elemID, localFace] = it->second;

			// input solver tag into rng
			mesh.data.rng[elemID*fpe + localFace] = solverTag;

		}

	}

}

std::vector<Index> residuum::mesh::exchange::gmsh::MeshConverter::reorderConnectivity(const Index* conn, residuum::mesh::exchange::gmsh::ElementType type){

	using ET = mesh::exchange::gmsh::ElementType;

	switch (type) {

		case ET::QuadP1:
			return {conn[0], conn[1], conn[3], conn[2]};


		case ET::HexP1:
			return {conn[2], conn[6], conn[3], conn[7], conn[1], conn[5], conn[0], conn[4]};
		
		// TODO: impelment remaining element types
		case ET::TriP1:
		case ET::TetP1:
		case ET::TriP2:
		case ET::QuadP2:
		case ET::TetP2:
		case ET::HexP2:
			return std::vector<Index>(conn, conn + nodesPerElement(type));

		default:
			return std::vector<Index>(conn, conn + nodesPerElement(type));
	
	}

}
