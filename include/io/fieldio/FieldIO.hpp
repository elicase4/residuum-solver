#ifndef RESIDUUM_IO_FIELDIO_FIELDIO_HPP
#define RESIDUUM_IO_FIELDIO_FIELDIO_HPP

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/Types.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "io/visualization/VTKWriter.hpp"
#include "io/visualization/VTUWriter.hpp"
#include "io/utils/Binary.hpp"
#include "mesh/Mesh.hpp"
#include "topology/TopologicalDOF.hpp"

/* Binary Nodal Field Format Description (PNDF)

[HEADER]
  magic:    uint32  = 0x504E4446  ("PNDF")
  version:  uint32  = 1
  numNodes: uint64
  numDOFs:  uint32
  time:     Real
  meshTag:  uint64

[DATA]
  values: Real[numNodes * numDOFs]
*/

namespace residuum {
	namespace io {
		namespace fieldio {

			// one .pndf snapshot: node-major values plus the time recorded in its header
			struct NodalFieldSnapshot {
				std::vector<Real> values;
				Real time;
			}; // struct NodalFieldSnapshot

			class FieldIO {
			public:

				template<Index numDOFs>
				static void writeVTK(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, Real time, const Real* algField, const std::vector<std::string>& dofNames, const std::string& filename, visualization::VTKWriter::Format fmt = visualization::VTKWriter::Format::ASCII);

				template<Index numDOFs>
				static void writeVTU(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, Real time, const Real* algField, const std::vector<std::string>& dofNames, const std::vector<std::string>& dofUnits, const std::string& filename);

				template<Index numDOFs>
				static std::vector<Real> reconstructNodalField(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, Real time, const Real* algField);

				template<Index numDOFs>
				static void writeBinary(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, Real time, const Real* algField, const std::string& filename);

				template<Index numDOFs>
				static void writeBinaryRaw(const mesh::Mesh& mesh, Real time, const std::vector<Real>& nodalField, const std::string& filename);

				template<Index numDOFs>
				static NodalFieldSnapshot readBinary(const mesh::Mesh& mesh, const std::string& filename);

				static constexpr uint32_t PNDF_MAGIC = 0x504E4446u;
				static constexpr uint32_t PNDF_VERSION = 1u;

			private:

				// fingerprint over mesh identity
				static uint64_t computeMeshTag(const mesh::Mesh& mesh);

			}; // class FieldIO

		} // namespace fieldio
	} // namespace io
} // namespace residuum

#include "FieldIO.tpp"

#endif
