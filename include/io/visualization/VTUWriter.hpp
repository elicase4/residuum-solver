#ifndef RESIDUUM_IO_VISUALIZATION_VTUWRITER_HPP
#define RESIDUUM_IO_VISUALIZATION_VTUWRITER_HPP

#include <fstream>
#include <string>

#include "core/Types.hpp"

namespace residuum {
	namespace io {
		namespace visualization {

			class VTUWriter {
			public:

				VTUWriter(const std::string& filename, Index numNodes, Index numElems);

				~VTUWriter();

				VTUWriter(const VTUWriter&) = delete;
				VTUWriter& operator=(const VTUWriter&) = delete;

				VTUWriter(VTUWriter&&) = default;
				VTUWriter& operator=(VTUWriter&&) = default;

				// section order: writePoints, writeCells, then point-data: begin, write*, end

				void writePoints(const Real* xyz, Index spatialDim);

				void writeCells(const Index* ien, Index nodesPerElem, int vtkCellType);

				void beginPointData();

				void writeScalar(const std::string& name, const Real* data, const std::string& unit = "");

				void writeVector(const std::string& name, const Real* data, Index numComponents, const std::string& unit = "");

				void endPointData();

			private:

				void close();

				std::ofstream ofs_;
				Index numNodes_;
				Index numElems_;

				enum class State { Open, PointsWritten, CellsWritten, InPointData, Done };
				State state_ = State::Open;

			}; // class VTUWriter

		} // namespace visualization
	} // namespace io
} // namespace residuum

#endif
