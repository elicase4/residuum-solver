#ifndef RESIDUUM_IO_GMSHREADER_HPP
#define RESIDUUM_IO_GMSHREADER_HPP

#include <cstdint>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "core/Types.hpp"
#include "io/utils/Binary.hpp"
#include "mesh/exchange/gmsh/Utils.hpp"
#include "mesh/exchange/gmsh/IntermediateMesh.hpp"

namespace residuum {
	namespace io {
		class GmshReader {
			public:
				
				static void read(mesh::exchange::gmsh::IntermediateMesh& mesh, const std::string& filename);

			private:
				
				enum class Format {
					ASCII,
					Binary
				}; // enum class Format

				struct VersionInfo {
					double version;
					Format format;
				};

				// high level version readers
				static VersionInfo readMeshFormat(std::istream& is);

				static void readMSH4(std::istream& is, mesh::exchange::gmsh::IntermediateMesh& mesh, Format fmt);
				
				static void readMSH2(std::istream& is, mesh::exchange::gmsh::IntermediateMesh& mesh, Format fmt);

				// internal helpers
				static void skipToSection(std::istream& is, const std::string& tag);

				static void deduceDimensions(mesh::exchange::gmsh::IntermediateMesh& mesh);

				static void readPhysicalNames(std::istream& is, std::unordered_map<Int, std::string>& names);

				static std::unordered_map<mesh::exchange::gmsh::EntityKey, Int, mesh::exchange::gmsh::EntityKeyHash> readEntities(std::istream& is, Format fmt);

				static void readNodes(std::istream& is, mesh::exchange::gmsh::IntermediateMesh& mesh, std::unordered_map<Index, Index>& tagToIdx, Format fmt);
				
				static void readElements(std::istream& is, mesh::exchange::gmsh::IntermediateMesh& mesh, const std::unordered_map<Index, Index>& tagToIdx, const std::unordered_map<mesh::exchange::gmsh::EntityKey, Int, mesh::exchange::gmsh::EntityKeyHash>& entityPhys, Format fmt);

		}; // class GmshReader
	} // namespace io
} // namespace residuum

#endif
