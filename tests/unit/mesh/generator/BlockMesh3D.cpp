#include <cmath>
#include <stdexcept>
#include <gtest/gtest.h>

#include "mesh/generator/BlockMesh3D.hpp"
#include "io/MeshIO.hpp"
#include "io/visualization/VTKWriter.hpp"

using namespace residuum::mesh::generator;
using namespace residuum::mesh;
using namespace residuum;

// ===================================================
// Basic Construction Tests
// ===================================================

TEST(BlockMesh3D, BasicConstruction) {
	// 2x2x2 unit cube with linear elements
	BlockMesh3D gen(2, 2, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1);
	Mesh mesh = gen.generate();

	EXPECT_TRUE(mesh.isValid());
	EXPECT_EQ(mesh.data.numNodes, 27);
	EXPECT_EQ(mesh.data.numElements, 8);
	EXPECT_EQ(mesh.data.nodesPerElement, 8);
	EXPECT_EQ(mesh.data.facesPerElement, 6);
	EXPECT_EQ(mesh.data.parametricDim, 3);
	EXPECT_EQ(mesh.data.spatialDim, 3);
	EXPECT_EQ(mesh.data.elementFamily, ElementFamily::Hex);
	EXPECT_EQ(mesh.data.basisType, BasisType::Lagrange);
}

TEST(BlockMesh3D, SingleElement) {
	// 1x1x1 unit cube with linear elements
	BlockMesh3D gen(1, 1, 1, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1);
	Mesh mesh = gen.generate();

	EXPECT_TRUE(mesh.isValid());
	EXPECT_EQ(mesh.data.numNodes, 8);
	EXPECT_EQ(mesh.data.numElements, 1);
}

TEST(BlockMesh3D, NonUnitCube) {
	// non-unit domain [2, 5] x [-1, 3] x [0, 4] with linear elements
	BlockMesh3D gen(3, 2, 1, 2.0, 5.0, -1.0, 3.0, 0.0, 4.0, 1, 1, 1);
	Mesh mesh = gen.generate();

	EXPECT_TRUE(mesh.isValid());
	EXPECT_EQ(mesh.data.numNodes, 24);
	EXPECT_EQ(mesh.data.numElements, 6);
}

TEST(BlockMesh3D, RejectsNonPositiveParameters) {
	EXPECT_THROW(BlockMesh3D(0, 2, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1).generate(), std::invalid_argument);
	EXPECT_THROW(BlockMesh3D(2, 0, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1).generate(), std::invalid_argument);
	EXPECT_THROW(BlockMesh3D(2, 2, 0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1).generate(), std::invalid_argument);
	EXPECT_THROW(BlockMesh3D(2, 2, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0, 1, 1).generate(), std::invalid_argument);
	EXPECT_THROW(BlockMesh3D(2, 2, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 0, 1).generate(), std::invalid_argument);
	EXPECT_THROW(BlockMesh3D(2, 2, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 0).generate(), std::invalid_argument);
}

// ===================================================
// Node Coordinate Tests
// ===================================================

TEST(BlockMesh3D, NodeCoordinates_UnitCube) {
	// 2x2x2 unit cube with linear elements
	BlockMesh3D gen(2, 2, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1);
	Mesh mesh = gen.generate();

	// check corner nodes
	const Real* node0 = mesh.getNodeCoord(0);
	EXPECT_DOUBLE_EQ(node0[0], 0.0);
	EXPECT_DOUBLE_EQ(node0[1], 0.0);
	EXPECT_DOUBLE_EQ(node0[2], 0.0);

	const Real* node2 = mesh.getNodeCoord(2);
	EXPECT_DOUBLE_EQ(node2[0], 1.0);
	EXPECT_DOUBLE_EQ(node2[1], 0.0);
	EXPECT_DOUBLE_EQ(node2[2], 0.0);

	const Real* node6 = mesh.getNodeCoord(6);
	EXPECT_DOUBLE_EQ(node6[0], 0.0);
	EXPECT_DOUBLE_EQ(node6[1], 1.0);
	EXPECT_DOUBLE_EQ(node6[2], 0.0);

	const Real* node18 = mesh.getNodeCoord(18);
	EXPECT_DOUBLE_EQ(node18[0], 0.0);
	EXPECT_DOUBLE_EQ(node18[1], 0.0);
	EXPECT_DOUBLE_EQ(node18[2], 1.0);

	const Real* node26 = mesh.getNodeCoord(26);
	EXPECT_DOUBLE_EQ(node26[0], 1.0);
	EXPECT_DOUBLE_EQ(node26[1], 1.0);
	EXPECT_DOUBLE_EQ(node26[2], 1.0);

	// check center node
	const Real* node13 = mesh.getNodeCoord(13);
	EXPECT_DOUBLE_EQ(node13[0], 0.5);
	EXPECT_DOUBLE_EQ(node13[1], 0.5);
	EXPECT_DOUBLE_EQ(node13[2], 0.5);
}

TEST(BlockMesh3D, NodeCoordinates_NonUniform) {
	// non-unit domain [2, 5] x [-1, 3] x [0, 4] with linear elements
	BlockMesh3D gen(3, 2, 1, 2.0, 5.0, -1.0, 3.0, 0.0, 4.0, 1, 1, 1);
	Mesh mesh = gen.generate();

	Real dx = (5.0 - 2.0) / 3.0;
	Real dz = (4.0 - 0.0) / 1.0;

	// bottom-front-left corner
	const Real* node0 = mesh.getNodeCoord(0);
	EXPECT_DOUBLE_EQ(node0[0], 2.0);
	EXPECT_DOUBLE_EQ(node0[1], -1.0);
	EXPECT_DOUBLE_EQ(node0[2], 0.0);

	// bottom-front-right corner
	const Real* node3 = mesh.getNodeCoord(3);
	EXPECT_DOUBLE_EQ(node3[0], 5.0);
	EXPECT_DOUBLE_EQ(node3[1], -1.0);
	EXPECT_DOUBLE_EQ(node3[2], 0.0);

	// bottom-back-left corner
	const Real* node8 = mesh.getNodeCoord(8);
	EXPECT_DOUBLE_EQ(node8[0], 2.0);
	EXPECT_DOUBLE_EQ(node8[1], 3.0);
	EXPECT_DOUBLE_EQ(node8[2], 0.0);

	// top-front-left corner
	const Real* node12 = mesh.getNodeCoord(12);
	EXPECT_DOUBLE_EQ(node12[0], 2.0);
	EXPECT_DOUBLE_EQ(node12[1], -1.0);
	EXPECT_DOUBLE_EQ(node12[2], dz);

	// check spacing
	const Real* node1 = mesh.getNodeCoord(1);
	EXPECT_DOUBLE_EQ(node1[0], 2.0 + dx);
	EXPECT_DOUBLE_EQ(node1[1], -1.0);
	EXPECT_DOUBLE_EQ(node1[2], 0.0);
}

TEST(BlockMesh3D, NodeOrdering) {
	// verify coordinate row-major ordering
	BlockMesh3D gen(2, 2, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1);
	Mesh mesh = gen.generate();

	for (Index k = 0; k <= 2; ++k){
		for (Index j = 0; j <= 2; ++j){
			for (Index i = 0; i <= 2; ++i){
				Index nodeID = (k * 9) + (j * 3) + i;
				const Real* coord = mesh.getNodeCoord(nodeID);
				EXPECT_DOUBLE_EQ(coord[0], 0.5 * i);
				EXPECT_DOUBLE_EQ(coord[1], 0.5 * j);
				EXPECT_DOUBLE_EQ(coord[2], 0.5 * k);
			}
		}
	}
}

// ===================================================
// Element Connectivity Tests
// ===================================================

TEST(BlockMesh3D, ElementConnectivity_SingleElement) {
	// 1x1x1 unit cube with a single element
	BlockMesh3D gen(1, 1, 1, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1);
	Mesh mesh = gen.generate();

	const Index* elem0 = mesh.getElementNodes(0);
	EXPECT_EQ(elem0[0], 0);
	EXPECT_EQ(elem0[1], 1);
	EXPECT_EQ(elem0[2], 2);
	EXPECT_EQ(elem0[3], 3);
	EXPECT_EQ(elem0[4], 4);
	EXPECT_EQ(elem0[5], 5);
	EXPECT_EQ(elem0[6], 6);
	EXPECT_EQ(elem0[7], 7);
}

TEST(BlockMesh3D, ElementConnectivity_2x2x2Mesh) {
	// 2x2x2 unit cube with linear elements
	BlockMesh3D gen(2, 2, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1);
	Mesh mesh = gen.generate();

	// element (0,0,0)
	const Index* elem0 = mesh.getElementNodes(0);
	EXPECT_EQ(elem0[0], 0);
	EXPECT_EQ(elem0[1], 1);
	EXPECT_EQ(elem0[2], 3);
	EXPECT_EQ(elem0[3], 4);
	EXPECT_EQ(elem0[4], 9);
	EXPECT_EQ(elem0[5], 10);
	EXPECT_EQ(elem0[6], 12);
	EXPECT_EQ(elem0[7], 13);

	// element (1,0,0)
	const Index* elem1 = mesh.getElementNodes(1);
	EXPECT_EQ(elem1[0], 1);
	EXPECT_EQ(elem1[1], 2);
	EXPECT_EQ(elem1[2], 4);
	EXPECT_EQ(elem1[3], 5);
	EXPECT_EQ(elem1[4], 10);
	EXPECT_EQ(elem1[5], 11);
	EXPECT_EQ(elem1[6], 13);
	EXPECT_EQ(elem1[7], 14);

	// element (0,0,1)
	const Index* elem4 = mesh.getElementNodes(4);
	EXPECT_EQ(elem4[0], 9);
	EXPECT_EQ(elem4[1], 10);
	EXPECT_EQ(elem4[2], 12);
	EXPECT_EQ(elem4[3], 13);
	EXPECT_EQ(elem4[4], 18);
	EXPECT_EQ(elem4[5], 19);
	EXPECT_EQ(elem4[6], 21);
	EXPECT_EQ(elem4[7], 22);
}

TEST(BlockMesh3D, AllElementsHaveValidNodes) {
	BlockMesh3D gen(3, 3, 3, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1);
	Mesh mesh = gen.generate();

	// Check all elements have valid node indices
	for (Index e = 0; e < mesh.data.numElements; ++e) {
		const Index* nodes = mesh.getElementNodes(e);
		for (Index a = 0; a < mesh.data.nodesPerElement; ++a) {
			EXPECT_LT(nodes[a], mesh.data.numNodes)
				<< "Element " << e << " node " << a << " out of bounds";
		}
	}
}

// ==================================================================
// Boundary Tag Tests
// ==================================================================

TEST(BlockMesh3D, BoundaryTags_SingleElement) {
	BlockMesh3D gen(1, 1, 1, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1);
	Mesh mesh = gen.generate();

	// Single element - all faces on boundary
	// LEFT=0, RIGHT=1, FRONT=2, BACK=3, BOTTOM=4, TOP=5
	EXPECT_TRUE(mesh.isOnBoundaryTag(0, 0, 0));  // LEFT
	EXPECT_TRUE(mesh.isOnBoundaryTag(0, 1, 1));  // RIGHT
	EXPECT_TRUE(mesh.isOnBoundaryTag(0, 2, 2));  // FRONT
	EXPECT_TRUE(mesh.isOnBoundaryTag(0, 3, 3));  // BACK
	EXPECT_TRUE(mesh.isOnBoundaryTag(0, 4, 4));  // BOTTOM
	EXPECT_TRUE(mesh.isOnBoundaryTag(0, 5, 5));  // TOP
}

TEST(BlockMesh3D, BoundaryTags_2x2x2Mesh) {
	BlockMesh3D gen(2, 2, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1);
	Mesh mesh = gen.generate();

	// Element 0 (ele_x=0, ele_y=0, ele_z=0): LEFT, FRONT, BOTTOM on boundary
	EXPECT_TRUE(mesh.isOnBoundaryTag(0, 0, 0));   // LEFT
	EXPECT_FALSE(mesh.isOnBoundary(0, 1));        // RIGHT (interior)
	EXPECT_TRUE(mesh.isOnBoundaryTag(0, 2, 2));   // FRONT
	EXPECT_FALSE(mesh.isOnBoundary(0, 3));        // BACK (interior)
	EXPECT_TRUE(mesh.isOnBoundaryTag(0, 4, 4));   // BOTTOM
	EXPECT_FALSE(mesh.isOnBoundary(0, 5));        // TOP (interior)

	// Element 7 (ele_x=1, ele_y=1, ele_z=1): RIGHT, BACK, TOP on boundary
	EXPECT_FALSE(mesh.isOnBoundary(7, 0));        // LEFT (interior)
	EXPECT_TRUE(mesh.isOnBoundaryTag(7, 1, 1));   // RIGHT
	EXPECT_FALSE(mesh.isOnBoundary(7, 2));        // FRONT (interior)
	EXPECT_TRUE(mesh.isOnBoundaryTag(7, 3, 3));   // BACK
	EXPECT_FALSE(mesh.isOnBoundary(7, 4));        // BOTTOM (interior)
	EXPECT_TRUE(mesh.isOnBoundaryTag(7, 5, 5));   // TOP
}

TEST(BlockMesh3D, InteriorElementsNoExternalFaces) {
	// 3x3x3 mesh - check the fully interior element has no boundary faces
	BlockMesh3D gen(3, 3, 3, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1);
	Mesh mesh = gen.generate();

	Index interiorElem = 1 * (3 * 3) + 1 * 3 + 1;  // ele_z=1, ele_y=1, ele_x=1

	EXPECT_FALSE(mesh.isOnBoundary(interiorElem, 0));  // LEFT
	EXPECT_FALSE(mesh.isOnBoundary(interiorElem, 1));  // RIGHT
	EXPECT_FALSE(mesh.isOnBoundary(interiorElem, 2));  // FRONT
	EXPECT_FALSE(mesh.isOnBoundary(interiorElem, 3));  // BACK
	EXPECT_FALSE(mesh.isOnBoundary(interiorElem, 4));  // BOTTOM
	EXPECT_FALSE(mesh.isOnBoundary(interiorElem, 5));  // TOP
}

// ==================================================================
// Higher-Order Element Tests
// ==================================================================

TEST(BlockMesh3D, QuadraticElements) {
	// Quadratic elements: px=2, py=2, pz=2
	BlockMesh3D gen(2, 2, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2, 2, 2);
	Mesh mesh = gen.generate();

	EXPECT_TRUE(mesh.isValid());
	EXPECT_EQ(mesh.data.nodesPerElement, 27);
	EXPECT_EQ(mesh.data.numElements, 8);
	EXPECT_EQ(mesh.data.numNodes, 125);
	EXPECT_EQ(mesh.data.basisOrder[0], 2);
	EXPECT_EQ(mesh.data.basisOrder[1], 2);
	EXPECT_EQ(mesh.data.basisOrder[2], 2);
}

TEST(BlockMesh3D, AnisotropicOrder) {
	// Different orders per axis: px=2, py=1, pz=1
	BlockMesh3D gen(2, 2, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2, 1, 1);
	Mesh mesh = gen.generate();

	EXPECT_TRUE(mesh.isValid());
	EXPECT_EQ(mesh.data.nodesPerElement, 12);    // (2+1)*(1+1)*(1+1)
	EXPECT_EQ(mesh.data.numElements, 8);         // 2*2*2
	EXPECT_EQ(mesh.data.numNodes, 45);           // (2*2+1)*(2*1+1)*(2*1+1)
	EXPECT_EQ(mesh.data.basisOrder[0], 2);
	EXPECT_EQ(mesh.data.basisOrder[1], 1);
	EXPECT_EQ(mesh.data.basisOrder[2], 1);
}

TEST(BlockMesh3D, HigherOrderConnectivity) {
	// Verify connectivity for a single quadratic element
	BlockMesh3D gen(1, 1, 1, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2, 2, 2);
	Mesh mesh = gen.generate();

	EXPECT_EQ(mesh.data.numElements, 1);
	EXPECT_EQ(mesh.data.numNodes, 27);

	const Index* nodes = mesh.getElementNodes(0);
	// Verify all 27 nodes are unique and valid
	for (Index i = 0; i < 27; ++i) {
		EXPECT_LT(nodes[i], 27);
		for (Index j = i + 1; j < 27; ++j) {
			EXPECT_NE(nodes[i], nodes[j])
				<< "Duplicate node in element: " << nodes[i];
		}
	}
}

TEST(BlockMesh3D, HigherOrderNodeCoordinates) {
	// Single quadratic element on the unit cube: the node grid must be the full
	// 3x3x3 fine grid (spacing 0.5), mirroring the BlockMesh2D regression test
	BlockMesh3D gen(1, 1, 1, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2, 2, 2);
	Mesh mesh = gen.generate();

	ASSERT_EQ(mesh.data.numNodes, 27);

	for (Index k = 0; k <= 2; ++k) {
		for (Index j = 0; j <= 2; ++j) {
			for (Index i = 0; i <= 2; ++i) {
				Index nodeID = (k * 9) + (j * 3) + i;
				const Real* coord = mesh.getNodeCoord(nodeID);
				EXPECT_DOUBLE_EQ(coord[0], 0.5 * i) << "node " << nodeID;
				EXPECT_DOUBLE_EQ(coord[1], 0.5 * j) << "node " << nodeID;
				EXPECT_DOUBLE_EQ(coord[2], 0.5 * k) << "node " << nodeID;
			}
		}
	}
}

// ==================================================================
// Mesh Validation Tests
// ==================================================================

TEST(BlockMesh3D, ValidationAfterGeneration) {
	BlockMesh3D gen(3, 4, 2, -1.0, 2.0, 0.0, 5.0, 0.0, 3.0, 1, 1, 1);
	Mesh mesh = gen.generate();

	EXPECT_TRUE(mesh.isValid());
}

TEST(BlockMesh3D, ClearMesh) {
	BlockMesh3D gen(2, 2, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1);
	Mesh mesh = gen.generate();

	EXPECT_TRUE(mesh.isValid());

	mesh.clear();

	EXPECT_EQ(mesh.data.numNodes, 0);
	EXPECT_EQ(mesh.data.numElements, 0);
	EXPECT_EQ(mesh.data.xyz.size(), 0);
	EXPECT_EQ(mesh.data.ien.size(), 0);
	EXPECT_FALSE(mesh.isValid());
}
