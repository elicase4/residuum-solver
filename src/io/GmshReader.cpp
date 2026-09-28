#include "io/GmshReader.hpp"

void residuum::io::GmshReader::skipToSection(std::istream& is, const std::string& tag) {

	std::string line;
	while (std::getline(is,line)) {
		if (line == tag) return;
	}
	
	throw std::runtime_error("GmshReader: section not found: " + tag);

}

void residuum::io::GmshReader::deduceDimensions(residuum::mesh::exchange::gmsh::IntermediateMesh& mesh) {

	Index maxParam = 0;
	for (const auto& eb : mesh.elementBlocks) {
		maxParam = std::max(maxParam, mesh::exchange::gmsh::parametricDimension(eb.type));
	}
	mesh.parametricDim = maxParam;

	const std::size_t numNodes = mesh.xyz.size() / 3;
	
	bool allZZero = true;
	
	for (std::size_t n = 0; n < numNodes && allZZero; ++n) {
		if (std::abs(mesh.xyz[n*3 + 2]) > 1e-14) {
			allZZero = false;
		}
	}
	mesh.spatialDim = allZZero ? 2 : 3;

}

void residuum::io::GmshReader::readPhysicalNames(std::istream& is, std::unordered_map<Int, std::string>& names) {

	int numGroups;
	is >> numGroups;
	
	for (int i = 0; i < numGroups; ++i) {
		
		int dim, tag;
		is >> dim >> tag;
		is >> std::ws;

		std::string name;
		std::getline(is, name);
		
		if (!name.empty() && name.front() == '"') {
			name = name.erase(0,1);
		}

		if (!name.empty() && name.back() == '"') {
			name.pop_back();
		}

		names[tag] = name;
	
	}

}

std::unordered_map<residuum::mesh::exchange::gmsh::EntityKey, Int, residuum::mesh::exchange::gmsh::EntityKeyHash> residuum::io::GmshReader::readEntities(std::istream& is, residuum::io::GmshReader::Format fmt) {

	std::unordered_map<residuum::mesh::exchange::gmsh::EntityKey, Int, residuum::mesh::exchange::gmsh::EntityKeyHash> entityPhys;

	std::size_t nPoints, nCurves, nSurfaces, nVolumes;

	if (fmt == residuum::io::GmshReader::Format::ASCII) {
		is >> nPoints >> nCurves >> nSurfaces >> nVolumes;
	} else {
		nPoints = residuum::io::binary::readLE<int64_t>(is);
		nCurves = residuum::io::binary::readLE<int64_t>(is);
		nSurfaces = residuum::io::binary::readLE<int64_t>(is);
		nVolumes = residuum::io::binary::readLE<int64_t>(is);
	}

	auto readEntryASCII = [&](Int entityDim, bool hasBox) {

		int tag;
		is >> tag;
		double a, b, c;
		is >> a >> b >> c;
		if (hasBox) {
			double d, e, f;
			is >> d >> e >> f;
		}

		std::size_t numPhys;
		is >> numPhys;
		for (std::size_t i = 0; i < numPhys; ++i) {
			int ptag;
			is >> ptag;
			if (i == 0) entityPhys[{entityDim, tag}] = ptag;
		}

		if (hasBox) {
			std::size_t numBound;
			is >> numBound;
			for (std::size_t i = 0; i < numBound; ++i){
				int btag;
				is >> btag;
			}
		}

	}; // readEntryASCII

	auto readEntryBinary = [&](Int entityDim, bool hasBox) {

		int tag = residuum::io::binary::readLE<int32_t>(is);
		double a = residuum::io::binary::readLE<double>(is);
		double b = residuum::io::binary::readLE<double>(is);
		double c = residuum::io::binary::readLE<double>(is);
		(void) a; (void) b; (void) c;

		if (hasBox) {
			residuum::io::binary::readLE<double>(is); residuum::io::binary::readLE<double>(is); residuum::io::binary::readLE<double>(is);
		}

		auto numPhys = residuum::io::binary::readLE<int64_t>(is);
		for (int64_t i = 0; i < numPhys; ++i) {
			int ptag = residuum::io::binary::readLE<int32_t>(is);
			if (i == 0) entityPhys[{entityDim, tag}] = ptag;
		}

		if (hasBox) {
			auto numBound = residuum::io::binary::readLE<int64_t>(is);
			for (int64_t i = 0; i < numBound; ++i) {
				residuum::io::binary::readLE<int32_t>(is);
			}
		}

	}; // readEntryBinary

	if (fmt == residuum::io::GmshReader::Format::ASCII) {
		for (std::size_t i = 0; i < nPoints; ++i) readEntryASCII(0, false);
		for (std::size_t i = 0; i < nCurves; ++i) readEntryASCII(1, true);
		for (std::size_t i = 0; i < nSurfaces; ++i) readEntryASCII(2, true);
		for (std::size_t i = 0; i < nVolumes; ++i) readEntryASCII(3, true);
	} else {
		for (std::size_t i = 0; i < nPoints; ++i) readEntryBinary(0, false);
		for (std::size_t i = 0; i < nCurves; ++i) readEntryBinary(1, true);
		for (std::size_t i = 0; i < nSurfaces; ++i) readEntryBinary(2, true);
		for (std::size_t i = 0; i < nVolumes; ++i) readEntryBinary(3, true);
	}

	return entityPhys;

}

void residuum::io::GmshReader::readNodes(std::istream& is, residuum::mesh::exchange::gmsh::IntermediateMesh& mesh, std::unordered_map<Index, Index>& tagToIdx, residuum::io::GmshReader::Format fmt) {

	std::size_t numBlocks, numNodes, minTag, maxTag;
	if (fmt == residuum::io::GmshReader::Format::ASCII) {
		is >> numBlocks >> numNodes >> minTag >> maxTag;
	} else {
		numBlocks = residuum::io::binary::readLE<int64_t>(is);
		numNodes = residuum::io::binary::readLE<int64_t>(is);
		minTag = residuum::io::binary::readLE<int64_t>(is);
		maxTag = residuum::io::binary::readLE<int64_t>(is);
	}

	tagToIdx.reserve(numNodes);
	mesh.xyz.reserve(numNodes * 3);

	for (std::size_t b = 0; b < numBlocks; ++b) {
		
		int entityDim, entityTag, parametric;
		std::size_t blockNodes;
		if (fmt == residuum::io::GmshReader::Format::ASCII) {
			is >> entityDim >> entityTag >> parametric >> blockNodes;
		} else {
			entityDim = residuum::io::binary::readLE<int32_t>(is);
			entityTag = residuum::io::binary::readLE<int32_t>(is);
			parametric = residuum::io::binary::readLE<int32_t>(is);
			blockNodes = residuum::io::binary::readLE<int64_t>(is);
		}
	
		residuum::mesh::exchange::gmsh::NodeBlock nb;
		nb.entityDim = entityDim;
		nb.entityTag = entityTag;
		nb.nodeIDs.resize(blockNodes);

		for (std::size_t n = 0; n < blockNodes; ++n){
			
			Index tag;
			
			if (fmt == residuum::io::GmshReader::Format::ASCII) {
				is >> tag;
			} else {
				tag = residuum::io::binary::readLE<int64_t>(is);
			}

			nb.nodeIDs[n] = tag;
		
		}

		for (std::size_t n = 0; n < blockNodes; ++n){
			
			double x, y, z;
			
			if (fmt == residuum::io::GmshReader::Format::ASCII) {
				is >> x >> y >> z;
			} else {
				x = residuum::io::binary::readLE<double>(is);
				y = residuum::io::binary::readLE<double>(is);
				z = residuum::io::binary::readLE<double>(is);
			}
			
			Index idx = mesh.xyz.size() / 3;
			mesh.xyz.push_back(static_cast<Real>(x));
			mesh.xyz.push_back(static_cast<Real>(y));
			mesh.xyz.push_back(static_cast<Real>(z));

			tagToIdx[nb.nodeIDs[n]] = idx;

		}

		for (auto& id : nb.nodeIDs) {
			id = tagToIdx[id];
		}

		mesh.nodeBlocks.push_back(std::move(nb));

	}

}

void residuum::io::GmshReader::readElements(std::istream& is, residuum::mesh::exchange::gmsh::IntermediateMesh& mesh, const std::unordered_map<Index, Index>& tagToIdx, const std::unordered_map<residuum::mesh::exchange::gmsh::EntityKey, Int, residuum::mesh::exchange::gmsh::EntityKeyHash>& entityPhys, residuum::io::GmshReader::Format fmt) {

	std::size_t numBlocks, numElements, minTag, maxTag;
	if (fmt == residuum::io::GmshReader::Format::ASCII) {
		is >> numBlocks >> numElements >> minTag >> maxTag;
	} else {
		numBlocks = residuum::io::binary::readLE<int64_t>(is);
		numElements = residuum::io::binary::readLE<int64_t>(is);
		minTag = residuum::io::binary::readLE<int64_t>(is);
		maxTag = residuum::io::binary::readLE<int64_t>(is);
	}

	for (std::size_t b = 0; b < numBlocks; ++b){

		int entityDim, entityTag, elemTypeInt;
		std::size_t blockElems;
		if (fmt == residuum::io::GmshReader::Format::ASCII) {
			is >> entityDim >> entityTag >> elemTypeInt >> blockElems;
		} else {
			entityDim = residuum::io::binary::readLE<int32_t>(is);
			entityTag = residuum::io::binary::readLE<int32_t>(is);
			elemTypeInt = residuum::io::binary::readLE<int32_t>(is);
			blockElems = residuum::io::binary::readLE<int64_t>(is);
		}

		residuum::mesh::exchange::gmsh::ElementBlock eb;
		eb.entityDim = entityDim;
		eb.entityTag = entityTag;
		eb.type = residuum::mesh::exchange::gmsh::elementTypeFromGmsh(elemTypeInt);
		eb.nodesPerElement = residuum::mesh::exchange::gmsh::nodesPerElement(eb.type);

		auto it = entityPhys.find({entityDim, entityTag});
		eb.physicalTag = (it != entityPhys.end()) ? (it->second) : (-1);

		eb.elementIDs.reserve(blockElems);
		eb.connectivity.reserve(blockElems * eb.nodesPerElement);

		for (std::size_t e = 0; e < blockElems; ++e) {

			Index elemTag;
			if (fmt == residuum::io::GmshReader::Format::ASCII) {
				is >> elemTag;
			} else {
				elemTag = residuum::io::binary::readLE<int64_t>(is);
			}
			eb.elementIDs.push_back(elemTag);

			std::vector<Index> raw(eb.nodesPerElement);
			for (Index n = 0; n < eb.nodesPerElement; ++n){
				Index nodeTag;
				if (fmt == residuum::io::GmshReader::Format::ASCII) {
					is >> nodeTag;
				} else {
					nodeTag = residuum::io::binary::readLE<int64_t>(is);
				}
				raw[n] = tagToIdx.at(nodeTag);
			}

			for (Index idx : raw){
				eb.connectivity.push_back(idx);
			}

		}

		mesh.elementBlocks.push_back(std::move(eb));

	}

}

residuum::io::GmshReader::VersionInfo residuum::io::GmshReader::readMeshFormat(std::istream& is) {

	double version;
	int fileType; // 0 = ASCII, 1 = binary
	int dataSize;
	is >> version >> fileType >> dataSize;

	residuum::io::GmshReader::VersionInfo vi;
	vi.version = version;
	vi.format = ((fileType == 0) ? (residuum::io::GmshReader::Format::ASCII) : (residuum::io::GmshReader::Format::Binary));

	if (vi.format == residuum::io::GmshReader::Format::Binary) {

		std::string rest;
		std::getline(is, rest);
		int32_t probe;
		is.read(reinterpret_cast<char*>(&probe), sizeof(probe));
		if (probe != 1){
			throw std::runtime_error("GmshReader: binary endian probe mismatch");
		}

	}

	return vi;

}

void residuum::io::GmshReader::read(residuum::mesh::exchange::gmsh::IntermediateMesh& mesh, const std::string& filename) {

	std::ifstream file(filename, std::ios::binary);
	if (!file.is_open()){
		throw std::runtime_error("GmshReader: cannot open file: " + filename);
	}

	skipToSection(file, "$MeshFormat");
	auto vi = readMeshFormat(file);

	if (vi.version >= 4.0){
		readMSH4(file, mesh, vi.format);
	} else {
		throw std::runtime_error("GmshReader: only MSH4 format supported.");
	}

}

void residuum::io::GmshReader::readMSH4(std::istream& is, residuum::mesh::exchange::gmsh::IntermediateMesh& mesh, residuum::io::GmshReader::Format fmt) {

	std::unordered_map<Index, Index> tagToIdx;
	std::unordered_map<residuum::mesh::exchange::gmsh::EntityKey, Int, residuum::mesh::exchange::gmsh::EntityKeyHash> entityPhys;

	std::string line;
	while(std::getline(is, line)) {
		if (line == "$PhysicalNames"){
			readPhysicalNames(is, mesh.physicalNames);
			std::getline(is, line);
		} else if (line == "$Entities") {
			entityPhys = readEntities(is, fmt);
			std::getline(is, line);
		} else if (line == "$Nodes") {
			readNodes(is, mesh, tagToIdx, fmt);
			std::getline(is, line);
		} else if (line == "$Elements") {
			readElements(is, mesh, tagToIdx, entityPhys, fmt);
			std::getline(is, line);
		}
	}

	deduceDimensions(mesh);

}
void residuum::io::GmshReader::readMSH2(std::istream&, residuum::mesh::exchange::gmsh::IntermediateMesh&, residuum::io::GmshReader::Format) {

	throw std::runtime_error("GmshReader: MSH2 not implemented");

}
