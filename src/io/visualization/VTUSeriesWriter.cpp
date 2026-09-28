#include "io/visualization/VTUSeriesWriter.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

residuum::io::visualization::VTUSeriesWriter::VTUSeriesWriter(std::string directory, std::string prefix, Index digitWidth) : directory_(std::move(directory)), prefix_(std::move(prefix)), digitWidth_(digitWidth) {}

std::string residuum::io::visualization::VTUSeriesWriter::addStep(Index step, Real time) {

	std::ostringstream stepStr;
	stepStr << std::setw(static_cast<int>(digitWidth_)) << std::setfill('0') << step;

	const std::string filename = prefix_ + "_" + stepStr.str() + ".vtu";
	entries_.emplace_back(filename, time);

	writePVD();

	return directory_ + "/" + filename;

}

void residuum::io::visualization::VTUSeriesWriter::writePVD() const {

	const std::string pvdPath = directory_ + "/" + prefix_ + ".pvd";

	std::ofstream ofs(pvdPath);
	if (!ofs.is_open()) {
		throw std::runtime_error("VTUSeriesWriter::writePVD: cannot open '" + pvdPath + "'");
	}

	ofs << "<?xml version=\"1.0\"?>\n";
	ofs << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
	ofs << "  <Collection>\n";

	for (const auto& [filename, time] : entries_) {
		ofs << "    <DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\"" << filename << "\"/>\n";
	}

	ofs << "  </Collection>\n";
	ofs << "</VTKFile>\n";

}
