#ifndef RESIDUUM_MESH_EXCHANGE_GMSH_UTILS_HPP
#define RESIDUUM_MESH_EXCHANGE_GMSH_UTILS_HPP

#include <vector>

#include "core/Types.hpp"
#include "mesh/BasisType.hpp"
#include "mesh/ElementFamily.hpp"
#include "mesh/exchange/gmsh/ElementType.hpp"

namespace residuum {
	namespace mesh {
		namespace exchange {
			namespace gmsh {

				// struct to handle non-unqiue tags between entity types
				struct EntityKey {
					int dim;
					int tag;
					
					// dimensions and tag must be equal
					bool operator==(const EntityKey& other) const {
						return ((dim == other.dim) && (tag == other.tag));
					}
				}; // struct EntityKey

				// hash for entities
				struct EntityKeyHash {
					std::size_t operator()(const EntityKey& k) const {
						return (std::hash<Int>{}(k.dim)^(std::hash<Int>{}(k.tag) << 1));
					}
				}; // struct EntityKeyHash

				// hash for index vectors
				struct VecHash {
					
					std::size_t operator()(const std::vector<Index>& v) const noexcept {
						
						std::size_t seed = v.size();
						for (auto x : v){
							seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
						}
						
						return seed;
					}

				}; // struct VecHash

				mesh::exchange::gmsh::ElementType elementTypeFromGmsh(int type);

				Index nodesPerElement(mesh::exchange::gmsh::ElementType type);

				Index parametricDimension(mesh::exchange::gmsh::ElementType type);

				Index facesPerElement(mesh::exchange::gmsh::ElementType type);

				std::vector<Index> basisOrder(mesh::exchange::gmsh::ElementType type);

				mesh::ElementFamily elementFamily(mesh::exchange::gmsh::ElementType type);

				mesh::BasisType basisType(mesh::exchange::gmsh::ElementType type);

				std::vector<Index> localFaceNodes(const Index* elemNodes, mesh::exchange::gmsh::ElementType type, Index face);

			} // namespace gmsh
		} // namespace exchange
	} // namespace mesh
} // namespace residuum

#endif
