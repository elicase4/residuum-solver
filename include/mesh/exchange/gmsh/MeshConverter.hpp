#ifndef RESIDUUM_MESH_EXCHANGE_GMSH_MESHCONVERTER_HPP
#define RESIDUUM_MESH_EXCHANGE_GMSH_MESHCONVERTER_HPP

#include <algorithm>
#include <stdexcept>
#include <unordered_map>

#include "core/Types.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/exchange/gmsh/IntermediateMesh.hpp"
#include "mesh/exchange/gmsh/Utils.hpp"

namespace residuum {
	namespace mesh {
		namespace exchange {
			namespace gmsh {

				class MeshConverter {
				public:
					
					static void toSolverMesh(mesh::Mesh& mesh, const IntermediateMesh& input, const std::unordered_map<Int, Int>& physicalGroupMap = {});

				private:

					static ElementType buildConnectivity(mesh::Mesh& mesh, const IntermediateMesh& input, std::unordered_map<Index, Index>& gmshToSolver);

					static void buildBoundaryTags(mesh::Mesh& mesh, const IntermediateMesh& input, ElementType cellType, std::unordered_map<Index, Index>& gmshToSolver, const std::unordered_map<Int, Int>& physicalGroupMap);

					static std::vector<Index> reorderConnectivity(const Index* conn, ElementType type);

				}; // class MeshConverter

			} // namespace gmsh
		} // namespace exchange
	} // namespace mesh
} // namespace residuum

#endif
