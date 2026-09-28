#include "application/mesh/MeshDispatcher.hpp"
#include "application/mesh/MeshConfig.hpp"

#include "application/mesh/GmshMeshGenerator.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/generator/BlockMesh2D.hpp"
#include "mesh/generator/BlockMesh3D.hpp"
#include "io/MeshIO.hpp"

#include "solver/logging/LoggerFactory.hpp"

#include <stdexcept>

namespace residuum {
	namespace application {
		namespace mesh {

			bool MeshDispatcher::run(const MeshConfig& config) {

				const auto logger = solver::logging::makeDriverLogger(config.logging.driver, "mesh");

				residuum::mesh::Mesh mesh;

				switch (config.type) {

					case MeshConfig::Type::Block2D: {

						const auto& b = config.block2D;

						residuum::mesh::generator::BlockMesh2D gen{b.nx, b.ny, b.xmin, b.xmax, b.ymin, b.ymax, b.Px, b.Py};
						mesh = gen.generate();

						logger.event("Block2D: " + std::to_string(b.nx) + "x" + std::to_string(b.ny) + " elements generated");
						break;
					}

					case MeshConfig::Type::Block3D: {

						const auto& b = config.block3D;

						residuum::mesh::generator::BlockMesh3D gen{b.nx, b.ny, b.nz, b.xmin, b.xmax, b.ymin, b.ymax, b.zmin, b.zmax, b.Px, b.Py, b.Pz};
						mesh = gen.generate();

						logger.event("Block3D: " + std::to_string(b.nx) + "x" + std::to_string(b.ny) + "x" + std::to_string(b.nz) + " elements generated");
						break;
					}

					case MeshConfig::Type::Gmsh: {

						if (config.inputFile.empty()) {
							throw std::runtime_error("MeshDispatcher: Gmsh type requires 'file' field in config");
						}

						// TODO: add physical group mapping argument
						residuum::application::mesh::GmshMeshGenerator gen{config.inputFile};
						mesh = gen.generate();

						logger.event("Gmsh import: " + config.inputFile + " -> " + std::to_string(mesh.data.numElements) + " elements");
						break;
					}

					default:
						throw std::runtime_error("MeshDispatcher: unknown mesh type");
				}

				if (!mesh.isValid()) {
					logger.error("generated mesh failed validity check");
					return false;
				}

				io::MeshIO::writeBinary(mesh, config.outputFile);

				logger.event("wrote " + config.outputFile);

				return true;
			}

		} // namespace mesh
	} // namespace application
} // namespace residuum
