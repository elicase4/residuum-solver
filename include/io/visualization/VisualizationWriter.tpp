#include "io/fieldio/FieldIO.hpp"

namespace residuum::io::visualization {

	template<Index numDOFs>
	std::string VisualizationWriter::writeField(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, Index step, Real time, const Real* algField, const std::vector<std::string>& dofNames, const std::vector<std::string>& dofUnits, const std::vector<fieldio::NamedNodalField>& extraFields) {

		// neither VTK nor a ParaView-rendered VTU surfaces a separate units attribute, so the unit is embedded directly in the array name, the only thing either viewer displays
		std::vector<std::string> combinedNames(dofNames.size());
		for (std::size_t c = 0; c < dofNames.size(); ++c) {
			combinedNames[c] = dofUnits[c].empty() ? dofNames[c] : (dofNames[c] + "[" + dofUnits[c] + "]");
		}

		if (format_ == Format::VTU) {

			const std::string filename = series_.addStep(step, time);
			fieldio::FieldIO::writeVTU<numDOFs>(mesh, topoDOF, bcRegistry, time, algField, combinedNames, filename, extraFields);
			return filename;

		}

		const std::string filename = directory_ + "/" + prefix_ + ".vtk";

		fieldio::FieldIO::writeVTK<numDOFs>(mesh, topoDOF, bcRegistry, time, algField, combinedNames, filename, visualization::VTKWriter::Format::ASCII, extraFields);
		return filename;

	}

} // namespace residuum::io::visualization
