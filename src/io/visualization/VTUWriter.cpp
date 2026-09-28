#include "io/visualization/VTUWriter.hpp"

#include <cassert>
#include <stdexcept>

residuum::io::visualization::VTUWriter::VTUWriter(const std::string& filename, Index numNodes, Index numElems) : numNodes_(numNodes), numElems_(numElems) {

	ofs_.open(filename, std::ios::out);
	if (!ofs_.is_open()) {
		throw std::runtime_error("VTUWriter: cannot open file '" + filename + "'");
	}

	ofs_ << "<?xml version=\"1.0\"?>\n";
	ofs_ << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	ofs_ << "  <UnstructuredGrid>\n";
	ofs_ << "    <Piece NumberOfPoints=\"" << numNodes_ << "\" NumberOfCells=\"" << numElems_ << "\">\n";

}

residuum::io::visualization::VTUWriter::~VTUWriter() {
	if (ofs_.is_open()) close();
}

void residuum::io::visualization::VTUWriter::writePoints(const Real* xyz, Index spatialDim) {

	assert(state_ == State::Open);

	if (spatialDim != 2 && spatialDim != 3) {
		throw std::runtime_error("VTUWriter::writePoints: spatialDim must be 2 or 3");
	}

	ofs_ << "      <Points>\n";
	ofs_ << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

	for (Index n = 0; n < numNodes_; ++n) {
		const Real* p = xyz + n * spatialDim;
		ofs_ << "          " << p[0] << ' ' << p[1] << ' ' << (spatialDim == 3 ? p[2] : Real(0)) << '\n';
	}

	ofs_ << "        </DataArray>\n";
	ofs_ << "      </Points>\n";

	state_ = State::PointsWritten;

}

void residuum::io::visualization::VTUWriter::writeCells(const Index* ien, Index nodesPerElem, int vtkCellType) {

	assert(state_ == State::PointsWritten);

	ofs_ << "      <Cells>\n";

	ofs_ << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
	for (Index e = 0; e < numElems_; ++e) {
		ofs_ << "          ";
		for (Index k = 0; k < nodesPerElem; ++k) {
			if (k) ofs_ << ' ';
			ofs_ << ien[e * nodesPerElem + k];
		}
		ofs_ << '\n';
	}
	ofs_ << "        </DataArray>\n";

	ofs_ << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
	ofs_ << "          ";
	for (Index e = 0; e < numElems_; ++e) {
		if (e) ofs_ << ' ';
		ofs_ << (e + 1) * nodesPerElem;
	}
	ofs_ << '\n';
	ofs_ << "        </DataArray>\n";

	ofs_ << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	ofs_ << "          ";
	for (Index e = 0; e < numElems_; ++e) {
		if (e) ofs_ << ' ';
		ofs_ << vtkCellType;
	}
	ofs_ << '\n';
	ofs_ << "        </DataArray>\n";

	ofs_ << "      </Cells>\n";

	state_ = State::CellsWritten;

}

void residuum::io::visualization::VTUWriter::beginPointData() {

	assert(state_ == State::CellsWritten);
	ofs_ << "      <PointData>\n";
	state_ = State::InPointData;

}

void residuum::io::visualization::VTUWriter::writeScalar(const std::string& name, const Real* data, const std::string& unit) {

	assert(state_ == State::InPointData);

	ofs_ << "        <DataArray type=\"Float64\" Name=\"" << name << "\"";
	if (!unit.empty()) ofs_ << " units=\"" << unit << "\"";
	ofs_ << " NumberOfComponents=\"1\" format=\"ascii\">\n";
	ofs_ << "          ";
	for (Index n = 0; n < numNodes_; ++n) {
		if (n) ofs_ << ' ';
		ofs_ << data[n];
	}
	ofs_ << '\n';
	ofs_ << "        </DataArray>\n";

}

void residuum::io::visualization::VTUWriter::writeVector(const std::string& name, const Real* data, Index numComponents, const std::string& unit) {

	assert(state_ == State::InPointData);

	if (numComponents < 2 || numComponents > 3) {
		throw std::runtime_error("VTUWriter::writeVector: numComponents must be 2 or 3");
	}

	ofs_ << "        <DataArray type=\"Float64\" Name=\"" << name << "\"";
	if (!unit.empty()) ofs_ << " units=\"" << unit << "\"";
	ofs_ << " NumberOfComponents=\"3\" format=\"ascii\">\n";

	for (Index n = 0; n < numNodes_; ++n) {
		ofs_ << "          ";
		for (Index c = 0; c < numComponents; ++c) {
			if (c) ofs_ << ' ';
			ofs_ << data[n * numComponents + c];
		}
		if (numComponents == 2) ofs_ << " 0.0";
		ofs_ << '\n';
	}

	ofs_ << "        </DataArray>\n";

}

void residuum::io::visualization::VTUWriter::endPointData() {

	assert(state_ == State::InPointData);
	ofs_ << "      </PointData>\n";
	state_ = State::CellsWritten; // allow more point/cell data sections if needed

}

void residuum::io::visualization::VTUWriter::close() {

	ofs_ << "    </Piece>\n";
	ofs_ << "  </UnstructuredGrid>\n";
	ofs_ << "</VTKFile>\n";
	ofs_.close();

	state_ = State::Done;

}
