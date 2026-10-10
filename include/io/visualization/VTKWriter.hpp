#ifndef RESIDUUM_IO_VISUALIZATION_VTKWRITER_HPP
#define RESIDUUM_IO_VISUALIZATION_VTKWRITER_HPP

#include <cassert>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/Types.hpp"
#include "io/utils/Binary.hpp"

namespace residuum {
	namespace io {
		namespace visualization {

			class VTKWriter {
			public:

				enum class Format { ASCII, Binary }; // enum class Format

				explicit VTKWriter(const std::string& filename, Format fmt = Format::ASCII);

				~VTKWriter();

				// remove copying
				VTKWriter(const VTKWriter&) = delete;
				VTKWriter& operator=(const VTKWriter&) = delete;

				// movable
				VTKWriter(VTKWriter&&) = default;
				VTKWriter& operator=(VTKWriter&&) = default;

				// writeHeader: 4 line VTK legacy file header (truncated to 256 chars)
				void writeHeader(const std::string& title = "solver output");

				// writePoints
				void writePoints(const Real* xyz, Index numNodes, Index spatialDim);

				// writeCells
				void writeCells(const Index* ien, Index numElems, Index nodesPerElem);

				// writeCellTypes: VTK cell type applies to each element (9 = VTK_QUAD - 2D quad4, 12 = VTK_HEXAHEDRON - 3D hex8)
				void writeCellTypes(int vtkType, Index numElems);

				// point-data sections order: begin, write data, end

				// beginPointData
				void beginPointData(Index numNodes);

				// writeScalar: data has length numNodes; name must embed its own unit (e.g. "Q[W/m^2]"), the only thing a viewer displays
				void writeScalar(const std::string& name, const Real* data, Index numNodes);

				// writeVector: data has length numNodes*numComponents, capped at 2 or 3 (a true tensor field needs a writeTensor); name must embed its unit, same as writeScalar
				void writeVector(const std::string& name, const Real* data, Index numNodes, Index numComponents);

				// endPointData
				void endPointData();

				// cell-data sections order: begin, write data, end

				// beginCellData
				void beginCellData(Index numElems);

				// writeScalarCell
				void writeScalarCell(const std::string& name, const Real* data, Index numElems);

				// endCellData
				void endCellData();

				// helpers

				// infer vtk cell types integer from spatialDim & nodesPerElement
				static int inferVTKCellType(Index spatialDim, Index nodesPerElement);

				// row-major to ccw ordering
				static std::vector<Index> rowMajorToCCW(const Index* rm, Index numNodes);

			private:

				std::ofstream ofs_;
				Format fmt_;

				// State to track ordering
				enum class State {
					Open, HeaderWritten, PointsWritten, CellsWritten,
					CellTypesWritten, InPointData, InCellData, Done
				}; // enum class State

				State state_ = State::Open;

				bool pointDataOpen_ = false;
				bool cellDataOpen_ = false;

			}; // class VTKWriter

		} // namespace visualization
	} // namespace io
} // namespace residuum

#endif
