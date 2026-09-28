#include "application/mesh/GmshMeshGenerator.hpp"

#include "io/GmshReader.hpp"
#include "mesh/exchange/gmsh/IntermediateMesh.hpp"
#include "mesh/exchange/gmsh/MeshConverter.hpp"

residuum::mesh::Mesh residuum::application::mesh::GmshMeshGenerator::generate(const std::unordered_map<Int, Int>& physicalGroupMap) const {

	residuum::mesh::exchange::gmsh::IntermediateMesh intermediateMesh;
	residuum::io::GmshReader::read(intermediateMesh, inputFile_);

	residuum::mesh::Mesh mesh;
	residuum::mesh::exchange::gmsh::MeshConverter::toSolverMesh(mesh, intermediateMesh, physicalGroupMap);

	return mesh;

}
