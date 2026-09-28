#include "io/visualization/VisualizationWriter.hpp"

#include <stdexcept>

residuum::io::visualization::VisualizationWriter::VisualizationWriter(std::string directory, std::string prefix, Format format, bool transient) : directory_(directory), prefix_(std::move(prefix)), format_(format), series_(std::move(directory), prefix_) {

	if (transient && format_ == Format::VTK) {
		throw std::runtime_error("VisualizationWriter: 'vtk' output format cannot produce a ParaView-readable time series (no .pvd support) -- use 'vtu' for a transient driver");
	}

}

residuum::utils::logging::CsvWriter residuum::io::visualization::VisualizationWriter::makeCsvWriter(const std::string& file, const std::vector<std::string>& columns) {

	return residuum::utils::logging::CsvWriter(file, columns);

}
