#ifndef RESIDUUM_IO_MESHIO_HPP
#define RESIDUUM_IO_MESHIO_HPP

#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/Types.hpp"
#include "io/visualization/VTKWriter.hpp"
#include "io/utils/Binary.hpp"
#include "mesh/Mesh.hpp"

/* Binary Mesh Format Description (PMSH)
[HEADER]
  magic:         uint32  = 0x504D5348  ("PMSH")
  version:       uint32  = 1
  parametricDim: uint32
  spatialDim:    uint32
  elementFamily: uint32  (mesh::ElementFamily)
  basisType:     uint32  (mesh::BasisType)
  numNodes:      uint64
  numElements:   uint64
  nodesPerElem:  uint32
  facesPerElem:  uint32
  numBasisOrders:uint32
  basisOrders:   uint32[numBasisOrders]
  hasExtractionOps: uint8  (0 = FEM, 1 = IGA)
  extractionOpSize: uint64 (0 if FEM)

[DATA]
  xyz:  Real[numNodes * spatialDim]       (floats/doubles, little-endian)
  ien:  Index[numElements * nodesPerElem] (uint64, little-endian)
  rng:  Int[numElements * facesPerElem]   (int32, little-endian)
  C:    Real[extractionOpSize]            (optional, IGA only)
*/

namespace residuum {
	namespace io {
		
		class MeshIO {
		public:

			// write mesh geometry to VTK legacy
			static void writeVTK(mesh::Mesh& mesh, const std::string& filename, visualization::VTKWriter::Format fmt = visualization::VTKWriter::Format::ASCII);

			// read/write binary PMSH
			static void writeBinary(const mesh::Mesh& mesh, const std::string& filename);
			static void readBinary(mesh::Mesh& mesh, const std::string& filename);

			static constexpr uint32_t PMSH_MAGIC = 0x504D5348u; // "PMSH" format identifier
			static constexpr uint32_t PMSH_VERSION = 1u;

		}; // class MeshIO
	
	} // namespace io
} // namespace residuum

#endif
