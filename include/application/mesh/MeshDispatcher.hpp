#ifndef RESIDUUM_APPLICATION_MESH_MESHDISPATCHER_HPP
#define RESIDUUM_APPLICATION_MESH_MESHDISPATCHER_HPP

#include "application/mesh/MeshConfig.hpp"

namespace residuum {
	namespace application {
		namespace mesh {

			class MeshDispatcher {
			public:

				static bool run(const MeshConfig& config);

			}; // class MeshDispatcher

		} // namespace mesh
	} // namespace application
} // namespace residuum

#endif
