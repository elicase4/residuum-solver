#include <filesystem>
#include <gtest/gtest.h>

#include "core/Types.hpp"
#include "io/GmshReader.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/exchange/gmsh/IntermediateMesh.hpp"

using namespace residuum;

TEST(GmshReader, QuadP1Structured) {

	mesh::exchange::gmsh::IntermediateMesh mesh;

	const std::filesystem::path input_mesh_path = std::filesystem::path(TEST_DATA_PATH) / "mesh/gmsh/msh/quad_p1.msh";

	io::GmshReader::read(mesh, input_mesh_path.string());

	// dimensions
	EXPECT_EQ(mesh.parametricDim, 2);
	EXPECT_EQ(mesh.spatialDim, 2);

	// node data
	EXPECT_EQ(mesh.xyz.size()/3, 9);
	EXPECT_EQ(mesh.nodeBlocks.size(), 9);

	// physical groups
	EXPECT_EQ(mesh.physicalNames.size(), 5);
	EXPECT_EQ(mesh.physicalNames.at(1), "bottom");
	EXPECT_EQ(mesh.physicalNames.at(2), "right");
	EXPECT_EQ(mesh.physicalNames.at(3), "top");
	EXPECT_EQ(mesh.physicalNames.at(4), "left");
	EXPECT_EQ(mesh.physicalNames.at(5), "domain");

	// element counts
	Index quadCount = 0;
	Index lineCount = 0;

	for (const auto& eb : mesh.elementBlocks) {

		if (eb.type == mesh::exchange::gmsh::ElementType::QuadP1) {
			quadCount += eb.elementIDs.size();

			EXPECT_EQ(eb.entityDim, 2);
			EXPECT_EQ(eb.physicalTag, 5);
			EXPECT_EQ(eb.nodesPerElement, 4);
		}

		if (eb.type == mesh::exchange::gmsh::ElementType::LineP1) {
			lineCount += eb.elementIDs.size();
			
			EXPECT_EQ(eb.entityDim, 1);
			EXPECT_GE(eb.physicalTag, 1);
			EXPECT_LE(eb.physicalTag, 4);
			EXPECT_EQ(eb.nodesPerElement, 2);
		}

	}

	EXPECT_EQ(quadCount, 4);
	EXPECT_EQ(lineCount, 8);

	// check coordinates
	EXPECT_DOUBLE_EQ(mesh.xyz[0], 0.0);
	EXPECT_DOUBLE_EQ(mesh.xyz[1], 0.0);
	EXPECT_DOUBLE_EQ(mesh.xyz[2], 0.0);

	// find mesh center in coordinate list
	bool foundCenter = false;

	for (Index n = 0; n < 9; ++n){
		
		const Real x = mesh.xyz[3*n + 0];
		const Real y = mesh.xyz[3*n + 1];

		if ((std::abs(x - 0.5) < 1e-13) && (std::abs(y - 0.5) < 1e-13)) {
			foundCenter = true;
		}

	}

	EXPECT_TRUE(foundCenter);

	// validate connectivity
	for (const auto& eb : mesh.elementBlocks){
		for (Index node : eb.connectivity) {
			EXPECT_LT(node, 9);
			EXPECT_GE(node, 0);
		}
	}

}

TEST(GmshReader, HexP1Structured) {

	mesh::exchange::gmsh::IntermediateMesh mesh;

	const std::filesystem::path input_mesh_path = std::filesystem::path(TEST_DATA_PATH) / "mesh/gmsh/msh/hex_p1.msh";

	io::GmshReader::read(mesh, input_mesh_path.string());

	// dimesnions
	EXPECT_EQ(mesh.parametricDim, 3);
	EXPECT_EQ(mesh.spatialDim, 3);

	// nodes
	EXPECT_EQ(mesh.xyz.size()/3, 27);
	EXPECT_EQ(mesh.nodeBlocks.size(), 27);

	// physical groups
	EXPECT_EQ(mesh.physicalNames.size(), 7);

	EXPECT_EQ(mesh.physicalNames.at(1), "left");
	EXPECT_EQ(mesh.physicalNames.at(2), "right");
	EXPECT_EQ(mesh.physicalNames.at(3), "front");
	EXPECT_EQ(mesh.physicalNames.at(4), "back");
	EXPECT_EQ(mesh.physicalNames.at(5), "bottom");
	EXPECT_EQ(mesh.physicalNames.at(6), "top");
	EXPECT_EQ(mesh.physicalNames.at(7), "domain");

	// elements
	Index hexCount = 0;
	Index quadFaceCount = 0;

	for (const auto& eb : mesh.elementBlocks) {

		if (eb.type == mesh::exchange::gmsh::ElementType::HexP1) {
			hexCount += eb.elementIDs.size();

			EXPECT_EQ(eb.entityDim, 3);
			EXPECT_EQ(eb.physicalTag, 7);
			EXPECT_EQ(eb.nodesPerElement, 8);
		}

		if (eb.type == mesh::exchange::gmsh::ElementType::QuadP1) {
			quadFaceCount += eb.elementIDs.size();

			EXPECT_EQ(eb.entityDim, 2);
			EXPECT_GE(eb.physicalTag, 1);
			EXPECT_LE(eb.physicalTag, 6);
			EXPECT_EQ(eb.nodesPerElement, 4);
		}

	}

	EXPECT_EQ(hexCount, 8);
	EXPECT_EQ(quadFaceCount, 24);

	// find mesh center in coordinate list
	bool foundCenter = false;

	for (Index n = 0; n < 27; ++n){
		
		const Real x = mesh.xyz[3*n + 0];
		const Real y = mesh.xyz[3*n + 1];
		const Real z = mesh.xyz[3*n + 2];

		if ((std::abs(x - 0.5) < 1e-13) && (std::abs(y - 0.5) < 1e-13) && (std::abs(z - 0.5)< 1e-13)) {
			foundCenter = true;
		}

	}

	EXPECT_TRUE(foundCenter);

	// check connectivity
	for (const auto& eb : mesh.elementBlocks){
		for (Index node : eb.connectivity){
			EXPECT_LT(node, 27);
			EXPECT_GE(node, 0);
		}
	}

}
