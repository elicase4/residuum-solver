#ifndef RESIDUUM_IO_VISUALIZATION_VISUALIZATIONWRITER_HPP
#define RESIDUUM_IO_VISUALIZATION_VISUALIZATIONWRITER_HPP

#include <string>
#include <vector>

#include "core/Types.hpp"

#include "fem/boundary/EssentialBoundaryRegistry.hpp"

#include "io/visualization/VTUSeriesWriter.hpp"

#include "mesh/Mesh.hpp"

#include "topology/TopologicalDOF.hpp"

#include "utils/logging/core/CsvWriter.hpp"

namespace residuum {
	namespace io {
		namespace visualization {

			class VisualizationWriter {
			public:

				enum class Format { VTK, VTU };

				VisualizationWriter(std::string directory, std::string prefix, Format format, bool transient);

				template<Index numDOFs>
				std::string writeField(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, Index step, Real time, const Real* algField, const std::vector<std::string>& dofNames, const std::vector<std::string>& dofUnits);

				static utils::logging::CsvWriter makeCsvWriter(const std::string& file, const std::vector<std::string>& columns);

			private:

				std::string directory_;
				std::string prefix_;
				Format format_;

				VTUSeriesWriter series_;

			}; // class VisualizationWriter

		} // namespace visualization
	} // namespace io
} // namespace residuum

#include "io/visualization/VisualizationWriter.tpp"

#endif
