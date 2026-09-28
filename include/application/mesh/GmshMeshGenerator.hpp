#ifndef RESIDUUM_APPLICATION_MESH_GMSHMESHGENERATOR_HPP
#define RESIDUUM_APPLICATION_MESH_GMSHMESHGENERATOR_HPP

#include <string>

#include "mesh/generator/MeshGenerator.hpp"

namespace residuum {
	namespace application {
		namespace mesh {

			class GmshMeshGenerator : public residuum::mesh::generator::MeshGenerator {
			public:

				explicit GmshMeshGenerator(std::string inputFile): inputFile_(std::move(inputFile)) {};

				residuum::mesh::Mesh generate(const std::unordered_map<Int, Int>& physicalGroupMap = {}) const override;

			private:

				std::string inputFile_;

			}; // class GmshMeshGenerator

		} // namespace mesh
	} // namespace application
} // namespace residuum

#endif
