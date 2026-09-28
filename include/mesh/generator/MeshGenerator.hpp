#ifndef RESIDUUM_MESH_GENERATOR_MESHGENERATOR_HPP
#define RESIDUUM_MESH_GENERATOR_MESHGENERATOR_HPP

#include <unordered_map>

#include "core/Types.hpp"
#include "mesh/Mesh.hpp"

namespace residuum {
	namespace mesh {
		namespace generator {

			class MeshGenerator {
			public:

				virtual ~MeshGenerator() = default;

				virtual Mesh generate(const std::unordered_map<Int, Int>& physicalGroupMap = {}) const = 0;

			}; // class MeshGenerator

		} // namespace generator
	} // namespace mesh
} // namespace residuum

#endif
