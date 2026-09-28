#ifndef RESIDUUM_MESH_EXCHANGE_GMSH_ELEMENTBLOCK_HPP
#define RESIDUUM_MESH_EXCHANGE_GMSH_ELEMENTBLOCK_HPP

#include <vector>

#include "core/Types.hpp"
#include "mesh/exchange/gmsh/ElementType.hpp"

namespace residuum {
	namespace mesh {
		namespace exchange {
			namespace gmsh {

				struct ElementBlock {

					Int entityDim = -1;
					Int entityTag = -1;

					ElementType type = ElementType::Unknown;

					Index nodesPerElement = 0;

					// flat connectivity
					std::vector<Index> connectivity;

					// file element IDs
					std::vector<Index> elementIDs;

					// physical group
					Int physicalTag = -1;

				}; // struct ElementBlock

			} // namespace gmsh
		} // namespace exchange
	} // namespace mesh
} // namespace residuum

#endif
