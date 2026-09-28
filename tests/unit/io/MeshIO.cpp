#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <iostream>

#include "core/Types.hpp"
#include "io/MeshIO.hpp"
#include "mesh/generator/BlockMesh2D.hpp"

using namespace residuum;
using namespace residuum::io;
using namespace residuum::mesh;

TEST(MeshIO, BinaryWriteRead_BlockMesh2D){

	const Index nx = 4;
	const Index ny = 4;

	generator::BlockMesh2D gen(nx, ny, 0.0, 1.0, 0.0, 1.0, 1, 1);
	Mesh mesh = gen.generate();

	ASSERT_TRUE(mesh.isValid());

	const std::filesystem::path output_path = std::filesystem::path(TEST_OUTPUT_PATH) / "test_mesh.pmsh";

	MeshIO::writeBinary(mesh, output_path.string());

	Mesh loaded;
	MeshIO::readBinary(loaded, output_path.string());

	EXPECT_TRUE(loaded.isValid());

	EXPECT_EQ(mesh.data.numNodes, loaded.data.numNodes);
	EXPECT_EQ(mesh.data.numElements, loaded.data.numElements);
	EXPECT_EQ(mesh.data.elementFamily, loaded.data.elementFamily);
	EXPECT_EQ(mesh.data.basisType, loaded.data.basisType);

	for (Index i = 0; i < mesh.data.xyz.size(); ++i){
		EXPECT_NEAR(mesh.data.xyz[i], loaded.data.xyz[i], 1e-14);
	}

	for (Index i = 0; i < mesh.data.ien.size(); ++i){
		EXPECT_NEAR(mesh.data.ien[i], loaded.data.ien[i], 1e-14);
	}

	for (Index i = 0; i < mesh.data.rng.size(); ++i){
		EXPECT_NEAR(mesh.data.rng[i], loaded.data.rng[i], 1e-14);
	}

}

TEST(MeshIO, WriteVTK_BlockMesh2D){

	const Index nx = 4;
	const Index ny = 4;

	mesh::generator::BlockMesh2D gen(nx, ny, 0.0, 1.0, 0.0, 1.0, 1, 1);
	Mesh mesh = gen.generate();

	ASSERT_TRUE(mesh.isValid());

	const std::filesystem::path output_path = std::filesystem::path(TEST_OUTPUT_PATH) / "test_mesh.vtk";
	MeshIO::writeVTK(mesh, output_path.string());

	std::ifstream file(output_path);
	EXPECT_TRUE(file.good());

	std::string line;

	std::getline(file, line);
	EXPECT_EQ(line, "# vtk DataFile Version 3.0");

	std::getline(file, line);
	EXPECT_EQ(line, "solver mesh");

	std::getline(file, line);
	EXPECT_EQ(line, "ASCII");

	std::getline(file, line);
	EXPECT_EQ(line, "DATASET UNSTRUCTURED_GRID");

	std::getline(file, line);
	EXPECT_EQ(line, "POINTS 25 double");

	file.seekg(0, std::ios::end);
	EXPECT_GT(file.tellg(), 0);

}
