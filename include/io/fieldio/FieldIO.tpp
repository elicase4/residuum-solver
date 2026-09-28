namespace residuum::io::fieldio {

	template<Index numDOFs>
	void FieldIO::writeVTK(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, Real time, const Real* algField, const std::vector<std::string>& dofNames, const std::string& filename, visualization::VTKWriter::Format fmt) {

		if (!mesh.isValid()){
			throw std::runtime_error("FieldIO::writeVTK: mesh is invalid");
		}

		if (dofNames.size() != topology::TopologicalDOF<numDOFs>::dofsPerNode){
			throw std::runtime_error("FieldIO::writeVTK: dofNames.size() (" + std::to_string(dofNames.size()) + ") must equal dofsPerNode (" + std::to_string(topology::TopologicalDOF<numDOFs>::dofsPerNode) + ")");
		}

		const int cellType = visualization::VTKWriter::inferVTKCellType(mesh.data.spatialDim, mesh.data.nodesPerElement);
		if (cellType == 0){
			throw std::runtime_error("FieldIO::writeVTK: unsupported spatialDim/nodesPerElement combination (" + std::to_string(mesh.data.spatialDim) + "D, " + std::to_string(mesh.data.nodesPerElement) + " nodes/elem");
		}

		// get ccw connectivity
		std::vector<Index> ienCCW(mesh.data.numElements * mesh.data.nodesPerElement);
		for (Index e = 0; e < mesh.data.numElements; ++e){

			std::vector<Index> ccw = visualization::VTKWriter::rowMajorToCCW(mesh.getElementNodes(e), mesh.data.nodesPerElement);

			for (Index k = 0; k < mesh.data.nodesPerElement; ++k){
				ienCCW[e * mesh.data.nodesPerElement + k] = ccw[k];
			}

		}

		// reconstruct full nodal field
		const std::vector<Real> nodalField = reconstructNodalField(mesh, topoDOF, bcRegistry, time, algField);

		// write file
		visualization::VTKWriter w(filename, fmt);
		w.writeHeader("solver field output");
		w.writePoints(mesh.data.xyz.data(), mesh.data.numNodes, mesh.data.spatialDim);
		w.writeCells(ienCCW.data(), mesh.data.numElements, mesh.data.nodesPerElement);
		w.writeCellTypes(cellType, mesh.data.numElements);

		// write point data
		w.beginPointData(mesh.data.numNodes);

		// extract interleaved data (always the case for topologicalDOF) from nodalField for stride-1 buffer to input to writeScalar
		std::vector<Real> componentBuf(mesh.data.numNodes);
		for (Index c = 0; c < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++c) {
			for (Index n = 0; n < mesh.data.numNodes; ++n){
				componentBuf[n] = nodalField[n*topology::TopologicalDOF<numDOFs>::dofsPerNode + c];
			}
			w.writeScalar(dofNames[c], componentBuf.data(), mesh.data.numNodes);
		}

		w.endPointData();

	}

	template<Index numDOFs>
	void FieldIO::writeVTU(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, Real time, const Real* algField, const std::vector<std::string>& dofNames, const std::vector<std::string>& dofUnits, const std::string& filename) {

		if (!mesh.isValid()){
			throw std::runtime_error("FieldIO::writeVTU: mesh is invalid");
		}

		if (dofNames.size() != topology::TopologicalDOF<numDOFs>::dofsPerNode){
			throw std::runtime_error("FieldIO::writeVTU: dofNames.size() (" + std::to_string(dofNames.size()) + ") must equal dofsPerNode (" + std::to_string(topology::TopologicalDOF<numDOFs>::dofsPerNode) + ")");
		}

		if (dofUnits.size() != dofNames.size()){
			throw std::runtime_error("FieldIO::writeVTU: dofUnits.size() (" + std::to_string(dofUnits.size()) + ") must equal dofNames.size() (" + std::to_string(dofNames.size()) + ")");
		}

		const int cellType = visualization::VTKWriter::inferVTKCellType(mesh.data.spatialDim, mesh.data.nodesPerElement);
		if (cellType == 0){
			throw std::runtime_error("FieldIO::writeVTU: unsupported spatialDim/nodesPerElement combination (" + std::to_string(mesh.data.spatialDim) + "D, " + std::to_string(mesh.data.nodesPerElement) + " nodes/elem");
		}

		// get ccw connectivity
		std::vector<Index> ienCCW(mesh.data.numElements * mesh.data.nodesPerElement);
		for (Index e = 0; e < mesh.data.numElements; ++e){

			std::vector<Index> ccw = visualization::VTKWriter::rowMajorToCCW(mesh.getElementNodes(e), mesh.data.nodesPerElement);

			for (Index k = 0; k < mesh.data.nodesPerElement; ++k){
				ienCCW[e * mesh.data.nodesPerElement + k] = ccw[k];
			}

		}

		// reconstruct full nodal field
		const std::vector<Real> nodalField = reconstructNodalField(mesh, topoDOF, bcRegistry, time, algField);

		// write file
		visualization::VTUWriter w(filename, mesh.data.numNodes, mesh.data.numElements);
		w.writePoints(mesh.data.xyz.data(), mesh.data.spatialDim);
		w.writeCells(ienCCW.data(), mesh.data.nodesPerElement, cellType);

		// write point data
		w.beginPointData();

		// extract interleaved data (always the case for topologicalDOF) from nodalField for stride-1 buffer to input to writeScalar
		std::vector<Real> componentBuf(mesh.data.numNodes);
		for (Index c = 0; c < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++c) {
			for (Index n = 0; n < mesh.data.numNodes; ++n){
				componentBuf[n] = nodalField[n*topology::TopologicalDOF<numDOFs>::dofsPerNode + c];
			}
			w.writeScalar(dofNames[c], componentBuf.data(), dofUnits[c]);
		}

		w.endPointData();

	}

	template<Index numDOFs>
	std::vector<Real> FieldIO::reconstructNodalField(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, Real time, const Real* algField) {

		// initialize containers for eval
		std::vector<Real> nodalField(mesh.data.numNodes * topology::TopologicalDOF<numDOFs>::dofsPerNode, Real(0.0));
		std::vector<Real> bcVal(topology::TopologicalDOF<numDOFs>::dofsPerNode, Real(0.0));

		// write free DOFs
		for (Index i = 0; i < topoDOF.numFreeDOFs(); ++i) {

			const Index topoDOFIdx = topoDOF.toTopological(i);
			const Index nodeId = topoDOF.getDOFNode(topoDOFIdx);
			const Index component = topoDOFIdx - (nodeId * topology::TopologicalDOF<numDOFs>::dofsPerNode);

			nodalField[nodeId * topology::TopologicalDOF<numDOFs>::dofsPerNode + component] = algField[i];

		}

		// write constrained DOFs
		for (Index topoDOFIdx = 0; topoDOFIdx < topoDOF.numGlobalDOFs(); ++topoDOFIdx) {

			if (!topoDOF.isConstrained(topoDOFIdx)) continue;

			const Index nodeId = topoDOF.getDOFNode(topoDOFIdx);
			const Index component = topoDOFIdx - (nodeId * topology::TopologicalDOF<numDOFs>::dofsPerNode);
			const Int tag = topoDOF.getConstraintTag(topoDOFIdx);
			const Real* xyz = mesh.getNodeCoord(nodeId);

			bool bcEntryFound = false;

			auto entries = bcRegistry.getEntries(tag);

			for (auto& entry : *entries) {
				entry->eval(time, xyz, bcVal.data());
				nodalField[nodeId * topology::TopologicalDOF<numDOFs>::dofsPerNode + component] = bcVal[component];
				bcEntryFound = true;
			}

			if (!bcEntryFound) {
				throw std::runtime_error("FieldIO::reconstructNodalField: failed to reconstruct constrained DOF.");
			}

		}

		return nodalField;

	}

	inline uint64_t FieldIO::computeMeshTag(const mesh::Mesh& mesh) {

		uint64_t hash = 14695981039346656037ull; // FNV-1a offset basis

		auto mix = [&hash](const void* data, std::size_t bytes) {
			const uint8_t* p = reinterpret_cast<const uint8_t*>(data);
			for (std::size_t i = 0; i < bytes; ++i) {
				hash ^= p[i];
				hash *= 1099511628211ull; // FNV-1a prime
			}
		};

		mix(mesh.data.xyz.data(), mesh.data.xyz.size() * sizeof(Real));
		mix(mesh.data.ien.data(), mesh.data.ien.size() * sizeof(Index));

		return hash;

	}

	template<Index numDOFs>
	void FieldIO::writeBinary(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, Real time, const Real* algField, const std::string& filename) {

		const std::vector<Real> nodalField = reconstructNodalField(mesh, topoDOF, bcRegistry, time, algField);
		FieldIO::writeBinaryRaw<numDOFs>(mesh, time, nodalField, filename);

	}

	template<Index numDOFs>
	void FieldIO::writeBinaryRaw(const mesh::Mesh& mesh, Real time, const std::vector<Real>& nodalField, const std::string& filename) {

		if (!mesh.isValid()) {
			throw std::runtime_error("FieldIO::writeBinaryRaw: mesh is invalid");
		}

		if (nodalField.size() != mesh.data.numNodes * numDOFs) {
			throw std::runtime_error("FieldIO::writeBinaryRaw: nodalField.size() (" + std::to_string(nodalField.size()) + ") must equal numNodes*numDOFs (" + std::to_string(mesh.data.numNodes * numDOFs) + ")");
		}

		std::ofstream ofs(filename, std::ios::binary);
		if (!ofs.is_open()) {
			throw std::runtime_error("FieldIO::writeBinaryRaw: cannot open file '" + filename + "'");
		}

		// header
		binary::writeLE<uint32_t>(ofs, PNDF_MAGIC);
		binary::writeLE<uint32_t>(ofs, PNDF_VERSION);
		binary::writeLE<uint64_t>(ofs, static_cast<uint64_t>(mesh.data.numNodes));
		binary::writeLE<uint32_t>(ofs, static_cast<uint32_t>(numDOFs));
		binary::writeLE<double>(ofs, static_cast<double>(time));
		binary::writeLE<uint64_t>(ofs, computeMeshTag(mesh));

		// data
		for (Real v : nodalField) {
			binary::writeLE<double>(ofs, static_cast<double>(v));
		}

		ofs.close();

	}

	template<Index numDOFs>
	NodalFieldSnapshot FieldIO::readBinary(const mesh::Mesh& mesh, const std::string& filename) {

		std::ifstream ifs(filename, std::ios::binary);
		if (!ifs.is_open()) {
			throw std::runtime_error("FieldIO::readBinary: cannot open file '" + filename + "'");
		}

		const uint32_t magic = binary::readLE<uint32_t>(ifs);
		if (magic != PNDF_MAGIC) {
			throw std::runtime_error("FieldIO::readBinary: bad magic number - file format is not a PNDF file");
		}

		const uint32_t version = binary::readLE<uint32_t>(ifs);
		if (version != PNDF_VERSION) {
			throw std::runtime_error("FieldIO::readBinary: unsupported PNDF version " + std::to_string(version));
		}

		const uint64_t numNodes = binary::readLE<uint64_t>(ifs);
		if (numNodes != static_cast<uint64_t>(mesh.data.numNodes)) {
			throw std::runtime_error("FieldIO::readBinary: numNodes mismatch (file has " + std::to_string(numNodes) + ", mesh has " + std::to_string(mesh.data.numNodes) + ")");
		}

		const uint32_t numDOFsFile = binary::readLE<uint32_t>(ifs);
		if (numDOFsFile != static_cast<uint32_t>(topology::TopologicalDOF<numDOFs>::dofsPerNode)) {
			throw std::runtime_error("FieldIO::readBinary: numDOFs mismatch (file has " + std::to_string(numDOFsFile) + ", expected " + std::to_string(topology::TopologicalDOF<numDOFs>::dofsPerNode) + ")");
		}

		const Real time = static_cast<Real>(binary::readLE<double>(ifs));

		const uint64_t meshTag = binary::readLE<uint64_t>(ifs);
		if (meshTag != computeMeshTag(mesh)) {
			throw std::runtime_error("FieldIO::readBinary: mesh tag mismatch - file '" + filename + "' was not written for the currently loaded mesh");
		}

		std::vector<Real> nodalField(numNodes * numDOFsFile);
		for (Real& v : nodalField) {
			v = static_cast<Real>(binary::readLE<double>(ifs));
		}

		return NodalFieldSnapshot{nodalField, time};

	}

} // namespace residuum::io::fieldio
