#ifndef RESIDUUM_APPLICATION_MESH_MESHAPPLICATION_HPP
#define RESIDUUM_APPLICATION_MESH_MESHAPPLICATION_HPP

#include "application/mesh/MeshConfig.hpp"

namespace residuum {
	namespace application {
		namespace mesh {

			class MeshApplication {
			public:

				explicit MeshApplication(const MeshConfig& config);
				
				int run();

			private:

				MeshConfig config_;

			}; // class MeshApplication

			// entry point called by the residuum subcommand dispatcher; argv[1] is the config path
			int runApp(int argc, char** argv);

		} // namespace mesh
	} // namespace application
} // namespace residuum

#endif
