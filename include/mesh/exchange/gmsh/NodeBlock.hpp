#ifndef RESIDUUM_MESH_EXCHANGE_GMSH_NODEBLOCK_HPP
#define RESIDUUM_MESH_EXCHANGE_GMSH_NODEBLOCK_HPP

#include <vector>

#include "core/Types.hpp"

namespace residuum {
	namespace mesh {
		namespace exchange {
			namespace gmsh {

				struct NodeBlock {

					Int entityDim = -1;
					Int entityTag = -1;

					std::vector<Index> nodeIDs;

				}; // struct NodeBlock

			} // namespace gmsh
		} // namespace exchange
	} // namespace mesh
} // namespace residuum

#endif
