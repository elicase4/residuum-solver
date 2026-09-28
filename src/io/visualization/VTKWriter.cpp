#include "io/visualization/VTKWriter.hpp"

residuum::io::visualization::VTKWriter::VTKWriter(const std::string& filename, Format fmt) : fmt_(fmt) {

	auto mode = std::ios::out;
	if (fmt_ == residuum::io::visualization::VTKWriter::Format::Binary){
		mode |= std::ios::binary;
	}
	ofs_.open(filename, mode);
	if (!ofs_.is_open()) {
		throw std::runtime_error("VTK Writer: cannot open file '" + filename + "'");
	}

}

residuum::io::visualization::VTKWriter::~VTKWriter() {
	if (ofs_.is_open()) ofs_.close();
}

void residuum::io::visualization::VTKWriter::writeHeader(const std::string& title) {

	assert(state_ == State::Open);

	ofs_ << "# vtk DataFile Version 3.0\n";
	ofs_ << title.substr(0, 256) << "\n";
	ofs_ << (fmt_ == Format::ASCII ? "ASCII\n" : "BINARY\n");
	ofs_ << "DATASET UNSTRUCTURED_GRID\n";

	state_ = State::HeaderWritten;

}

void residuum::io::visualization::VTKWriter::writePoints(const Real* xyz, Index numNodes, Index spatialDim) {

	assert(state_ == State::HeaderWritten);

	if (spatialDim != 2 && spatialDim != 3){
		throw std::runtime_error("VTKWriter::writePoints: spatialDim must be 2 or 3");
	}

	ofs_ << "POINTS " << numNodes << " double\n";

	for (Index n = 0; n < numNodes; ++n) {

		const double* p = xyz + n*spatialDim;

		if (fmt_ == Format::ASCII){
			ofs_ << p[0] << ' ' << p[1] << ' ';
			ofs_ << (spatialDim == 3 ? p[2] : static_cast<double>(0.0)) << '\n';
		} else {
			binary::writeBE<double>(ofs_, static_cast<double>(p[0]));
			binary::writeBE<double>(ofs_, static_cast<double>(p[1]));
			binary::writeBE<double>(ofs_, static_cast<double>( spatialDim == 3 ? p[2] : static_cast<Real>(0.0) ));
		}

	}

	if (fmt_ == Format::Binary) ofs_ << '\n';

	state_ = State::PointsWritten;

}

void residuum::io::visualization::VTKWriter::writeCells(const Index* ien, Index numElems, Index nodesPerElem) {

	assert(state_ == State::PointsWritten);

	const Index listSize = numElems * (nodesPerElem + 1);
	ofs_ << "CELLS " << numElems << ' ' << listSize << '\n';

	for (Index e = 0; e < numElems; ++e) {

		const Index* nodes = ien + e*nodesPerElem;

		if (fmt_ == Format::ASCII){
			ofs_ << nodesPerElem;
			for (Index k = 0; k < nodesPerElem; ++k){
				ofs_ << ' ' << nodes[k];
			}
			ofs_ << '\n';
		} else {
			binary::writeBE<uint32_t>(ofs_, static_cast<uint32_t>(nodesPerElem));
			for (Index k = 0; k < nodesPerElem; ++k){
				binary::writeBE<uint32_t>(ofs_, static_cast<uint32_t>(nodes[k]));
			}
		}

	}

	if (fmt_ == Format::Binary) ofs_ << '\n';

	state_ = State::CellsWritten;

}

void residuum::io::visualization::VTKWriter::writeCellTypes(int vtkType, Index numElems) {

	assert(state_ == State::CellsWritten);

	ofs_ << "CELL_TYPES " << numElems << '\n';

	for (Index e = 0; e < numElems; ++e) {

		if (fmt_ == Format::ASCII){
			ofs_ << vtkType << '\n';
		} else {
			binary::writeBE<uint32_t>(ofs_, static_cast<uint32_t>(vtkType));
		}

	}

	if (fmt_ == Format::Binary) ofs_ << '\n';

	state_ = State::CellTypesWritten;

}

void residuum::io::visualization::VTKWriter::beginPointData(Index numNodes) {

	assert(state_ == State::CellTypesWritten);
	ofs_ << "POINT_DATA " << numNodes << '\n';
	state_ = State::InPointData;
	pointDataOpen_ = true;

}

void residuum::io::visualization::VTKWriter::writeScalar(const std::string& name, const Real* data, Index numNodes) {

	assert(state_ == State::InPointData && pointDataOpen_);

	ofs_ << "SCALARS " << name << " double 1\n";
	ofs_ << "LOOKUP_TABLE default\n";

	for (Index n = 0; n < numNodes; ++n) {

		if (fmt_ == Format::ASCII) {
			ofs_ << data[n] << '\n';
		} else {
			binary::writeBE<double>(ofs_, static_cast<double>(data[n]));
		}

	}

	if (fmt_ == Format::Binary) ofs_ << '\n';

}

void residuum::io::visualization::VTKWriter::writeVector(const std::string& name, const Real* data, Index numNodes, Index numComponents) {

	assert(state_ == State::InPointData && pointDataOpen_);

	if (numComponents < 2 || numComponents > 3) {
		throw std::runtime_error("VTKWriter::writeVector: numComponents must be 2 or 3");
	}

	ofs_ << "VECTORS " << name << " double\n";

	for (Index n = 0; n < numNodes; ++n) {

		if (fmt_ == Format::ASCII) {
			for (Index c = 0; c < numComponents; ++c){
				if (c) ofs_ << ' ';
				ofs_ << data[n*numComponents + c];
			}
			if (numComponents == 2) ofs_ << " 0.0";
			ofs_ << '\n';
		} else {
			for (Index c = 0; c < numComponents; ++c){
				binary::writeBE<double>(ofs_, static_cast<double>(data[n*numComponents + c]));
			}
			if (numComponents == 2) binary::writeBE<double>(ofs_, 0.0);
		}

	}

	if (fmt_ == Format::Binary) ofs_ << '\n';

}

void residuum::io::visualization::VTKWriter::endPointData() {

	assert(pointDataOpen_);
	pointDataOpen_ = false;
	state_ = State::CellTypesWritten; // allow for more cell types to be written

}

void residuum::io::visualization::VTKWriter::beginCellData(Index numElems) {

	assert(state_ == State::CellTypesWritten && !cellDataOpen_);
	ofs_ << "CELL_DATA " << numElems << '\n';
	state_ = State::InCellData;
	cellDataOpen_ = true;

}

void residuum::io::visualization::VTKWriter::writeScalarCell(const std::string& name, const Real* data, Index numElems) {

	assert(state_ == State::InCellData && cellDataOpen_);
	ofs_ << "SCALARS " << name << " double 1\n";
	ofs_ << "LOOKUP_TABLE default\n";

	for (Index e = 0; e < numElems; ++e){

		if (fmt_ == Format::ASCII) {
			ofs_ << data[e] << '\n';
		} else {
			binary::writeBE<double>(ofs_, static_cast<double>(data[e]));
		}

	}

	if (fmt_ == Format::Binary) ofs_ << '\n';

}

void residuum::io::visualization::VTKWriter::endCellData() {

	assert(cellDataOpen_);
	cellDataOpen_ = false;
	state_ = State::Done;

}

int residuum::io::visualization::VTKWriter::inferVTKCellType(Index spatialDim, Index nodesPerElement) {

	if (spatialDim == 2) {
		switch (nodesPerElement) {
			case 3: return 5;  // VTK_TRIANGLE
			case 4: return 9;  // VTK_QUAD
			case 6: return 22; // VTK_QUADRATIC_TRIANGLE
			case 8: return 23; // VTK_QUADRATIC_QUAD
			default: return 0;
		}
	}

	if (spatialDim == 3) {
		switch (nodesPerElement) {
			case 4:  return 10; // VTK_TETRA
			case 8:  return 12; // VTK_HEXAHEDRON
			case 10: return 24; // VTK_QUADRATIC_TETRA
			case 20: return 25; // VTK_QUADRATIC_HEXAHEDRON
			default: return 0;
		}
	}

	return 0;
}

std::vector<Index> residuum::io::visualization::VTKWriter::rowMajorToCCW(const Index* rm, const Index numNodes) {

	switch (numNodes) {
		case 4:
			// 2D P1 quad
			return { rm[0], rm[1], rm[3], rm[2] };
		case 8:
			// 3D P1 quad
			return { rm[0], rm[1], rm[3], rm[2],
					 rm[4], rm[5], rm[7], rm[6] };
		default:
			throw std::runtime_error("VTKWriter::rowMajorToCCW: unsupported element size " + std::to_string(numNodes));
	}

}
