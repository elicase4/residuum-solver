#ifndef RESIDUUM_MESH_EXCHANGE_GMSH_INTERMEDIATEMESH_HPP
#define RESIDUUM_MESH_EXCHANGE_GMSH_INTERMEDIATEMESH_HPP

#include <string>
#include <unordered_map>
#include <vector>

#include "core/Types.hpp"

#include "mesh/exchange/gmsh/NodeBlock.hpp"
#include "mesh/exchange/gmsh/ElementBlock.hpp"

namespace residuum {
	namespace mesh {
		namespace exchange {
			namespace gmsh {

				struct IntermediateMesh {

					// dimensions
					Index parametricDim = 0;
					Index spatialDim = 0;

					// row-major node coordinates
					std::vector<Real> xyz;

					// mesh entities
					std::vector<NodeBlock> nodeBlocks;
					std::vector<ElementBlock> elementBlocks;

					// physical group metadata
					std::unordered_map<Int, std::string> physicalNames;
					
					// block accessors
					std::vector<const ElementBlock*> cellBlocks() const;
					std::vector<const ElementBlock*> boundaryBlocks() const;
					
					void clear() {
						parametricDim = 0; spatialDim = 0; xyz.clear(); nodeBlocks.clear(); elementBlocks.clear(); physicalNames.clear();
					}

					bool empty () const {
						return (xyz.empty() || elementBlocks.empty());
					}

				}; // struct IntermediateMesh

			} // namespace gmsh
		} // namespace exchange
	} // namespace mesh
} // namespace residuum

#endif
