#ifndef RESIDUUM_APPLICATION_MESH_MESHCONFIGPARSER_HPP
#define RESIDUUM_APPLICATION_MESH_MESHCONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include "application/mesh/MeshConfig.hpp"

namespace residuum {
	namespace application {
		namespace mesh {

			class MeshConfigParser {
			public:

				static MeshConfig::Type parseMeshType(const std::string& str);

				static MeshConfig read(const std::string& filename);
			
			}; // class MeshConfigReader

		} // namespace mesh
	} // namespace application
} // namespace residuum

#endif
