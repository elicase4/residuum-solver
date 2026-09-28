#include "io/fieldio/FieldIO.hpp"

namespace residuum::io::visualization {

	template<Index numDOFs>
	std::string VisualizationWriter::writeField(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, Index step, Real time, const Real* algField, const std::vector<std::string>& dofNames, const std::vector<std::string>& dofUnits) {

		if (format_ == Format::VTU) {

			const std::string filename = series_.addStep(step, time);
			fieldio::FieldIO::writeVTU<numDOFs>(mesh, topoDOF, bcRegistry, time, algField, dofNames, dofUnits, filename);
			return filename;

		}

		const std::string filename = directory_ + "/" + prefix_ + ".vtk";

		std::vector<std::string> combinedNames(dofNames.size());
		for (std::size_t c = 0; c < dofNames.size(); ++c) {
			combinedNames[c] = dofUnits[c].empty() ? dofNames[c] : (dofNames[c] + "[" + dofUnits[c] + "]");
		}

		fieldio::FieldIO::writeVTK<numDOFs>(mesh, topoDOF, bcRegistry, time, algField, combinedNames, filename);
		return filename;

	}

} // namespace residuum::io::visualization
