#include <cmath>
#include <stdexcept>
#include <gtest/gtest.h>

#include "mesh/generator/BlockMesh2D.hpp"
#include "io/MeshIO.hpp"
#include "io/visualization/VTKWriter.hpp"

using namespace residuum::mesh::generator;
using namespace residuum::mesh;
using namespace residuum;

// ===================================================
// Basic Construction Tests
// ===================================================

TEST(BlockMesh2D, BasicConstruction) {
	// 2x2 unit square with linear elements
	BlockMesh2D gen(2, 2, 0.0, 1.0, 0.0, 1.0, 1, 1);
	Mesh mesh = gen.generate();

	EXPECT_TRUE(mesh.isValid());
	EXPECT_EQ(mesh.data.numNodes, 9);
	EXPECT_EQ(mesh.data.numElements, 4);
	EXPECT_EQ(mesh.data.nodesPerElement, 4);
	EXPECT_EQ(mesh.data.facesPerElement, 4);
	EXPECT_EQ(mesh.data.parametricDim, 2);
	EXPECT_EQ(mesh.data.spatialDim, 2);
	EXPECT_EQ(mesh.data.elementFamily, ElementFamily::Quad);
	EXPECT_EQ(mesh.data.basisType, BasisType::Lagrange);
}

TEST(BlockMesh2D, SingleElement) {
	// 1x1 unit square with linear elements
	BlockMesh2D gen(1, 1, 0.0, 1.0, 0.0, 1.0, 1, 1);
	Mesh mesh = gen.generate();

	EXPECT_TRUE(mesh.isValid());
	EXPECT_EQ(mesh.data.numNodes, 4);
	EXPECT_EQ(mesh.data.numElements, 1);
}

TEST(BlockMesh2D, NonUnitSquare) {
	// non-unit domain [2, 5] x [-1, 3] with linear elements
	BlockMesh2D gen(3, 2, 2.0, 5.0, -1.0, 3.0, 1, 1);
	Mesh mesh = gen.generate();

	EXPECT_TRUE(mesh.isValid());
	EXPECT_EQ(mesh.data.numNodes, 12);
	EXPECT_EQ(mesh.data.numElements, 6);
}

TEST(BlockMesh2D, RejectsNonPositiveParameters) {
	EXPECT_THROW(BlockMesh2D(0, 2, 0.0, 1.0, 0.0, 1.0, 1, 1).generate(), std::invalid_argument);
	EXPECT_THROW(BlockMesh2D(2, 0, 0.0, 1.0, 0.0, 1.0, 1, 1).generate(), std::invalid_argument);
	EXPECT_THROW(BlockMesh2D(2, 2, 0.0, 1.0, 0.0, 1.0, 0, 1).generate(), std::invalid_argument);
	EXPECT_THROW(BlockMesh2D(2, 2, 0.0, 1.0, 0.0, 1.0, 1, 0).generate(), std::invalid_argument);
}

// ===================================================
// Node Coordinate Tests
// ===================================================

TEST(BlockMesh2D, NodeCoordinates_UnitSquare) {
	// 2x2 unit square with linear elements
	BlockMesh2D gen(2, 2, 0.0, 1.0, 0.0, 1.0, 1, 1);
	Mesh mesh = gen.generate();

	// check corner nodes
	const Real* node0 = mesh.getNodeCoord(0);
	EXPECT_DOUBLE_EQ(node0[0], 0.0);
	EXPECT_DOUBLE_EQ(node0[1], 0.0);

	const Real* node2 = mesh.getNodeCoord(2);
	EXPECT_DOUBLE_EQ(node2[0], 1.0);
	EXPECT_DOUBLE_EQ(node2[1], 0.0);

	const Real* node6 = mesh.getNodeCoord(6);
	EXPECT_DOUBLE_EQ(node6[0], 0.0);
	EXPECT_DOUBLE_EQ(node6[1], 1.0);

	const Real* node8 = mesh.getNodeCoord(8);
	EXPECT_DOUBLE_EQ(node8[0], 1.0);
	EXPECT_DOUBLE_EQ(node8[1], 1.0);

	// check center node
	const Real* node4 = mesh.getNodeCoord(4);
	EXPECT_DOUBLE_EQ(node4[0], 0.5);
	EXPECT_DOUBLE_EQ(node4[1], 0.5);
}

TEST(BlockMesh2D, NodeCoordinates_NonUniform) {
	// non-unit domain [2, 5] x [-1, 3] with linear elements
	BlockMesh2D gen(3, 2, 2.0, 5.0, -1.0, 3.0, 1, 1);
	Mesh mesh = gen.generate();

	Real dx = (5.0 - 2.0) / 3.0;
	Real dy = (3.0 - (-1.0)) / 2.0;

	// bottom-left corner
	const Real* node0 = mesh.getNodeCoord(0);
	EXPECT_DOUBLE_EQ(node0[0], 2.0);
	EXPECT_DOUBLE_EQ(node0[1], -1.0);

	// bottom-right corner
	const Real* node3 = mesh.getNodeCoord(3);
	EXPECT_DOUBLE_EQ(node3[0], 5.0);
	EXPECT_DOUBLE_EQ(node3[1], -1.0);

	// top-left corner
	const Real* node8 = mesh.getNodeCoord(8);
	EXPECT_DOUBLE_EQ(node8[0], 2.0);
	EXPECT_DOUBLE_EQ(node8[1], 3.0);

	// check spacing
	const Real* node1 = mesh.getNodeCoord(1);
	EXPECT_DOUBLE_EQ(node1[0], 2.0 + dx);
	EXPECT_DOUBLE_EQ(node1[1], -1.0);

	// check center node
	const Real* node4 = mesh.getNodeCoord(4);
	EXPECT_DOUBLE_EQ(node4[0], 2.0);
	EXPECT_DOUBLE_EQ(node4[1], -1.0 + dy);
}

TEST(BlockMesh2D, NodeOrdering) {
	// verify coordinate row-major ordering
	BlockMesh2D gen(2, 2, 0.0, 1.0, 0.0, 1.0, 1, 1);
	Mesh mesh = gen.generate();

	for (Index j = 0; j <= 2; ++j){
		for (Index i = 0; i <= 2; ++i){
			Index nodeID = (j * 3) + i;
			const Real* coord = mesh.getNodeCoord(nodeID);
			EXPECT_DOUBLE_EQ(coord[0], 0.5 * i);
			EXPECT_DOUBLE_EQ(coord[1], 0.5 * j);
		}
	}
}

// ===================================================
// Element Connectivity Tests
// ===================================================

TEST(BlockMesh2D, ElementConnectivity_SingleElement) {
	// 2x2 unit square with single elements
	BlockMesh2D gen(1, 1, 0.0, 1.0, 0.0, 1.0, 1, 1);
	Mesh mesh = gen.generate();

	// check connection of a single element
	const Index* elem0 = mesh.getElementNodes(0);
	EXPECT_EQ(elem0[0], 0);
	EXPECT_EQ(elem0[1], 1);
	EXPECT_EQ(elem0[2], 2);
	EXPECT_EQ(elem0[3], 3);
}

TEST(BlockMesh2D, ElementConnectivity_2x2Mesh) {
	// 2x2 unit square with linear elements
	BlockMesh2D gen(2, 2, 0.0, 1.0, 0.0, 1.0, 1, 1);
	Mesh mesh = gen.generate();

	// check connection of the elements
	const Index* elem0 = mesh.getElementNodes(0);
	EXPECT_EQ(elem0[0], 0);
	EXPECT_EQ(elem0[1], 1);
	EXPECT_EQ(elem0[2], 3);
	EXPECT_EQ(elem0[3], 4);

	const Index* elem1 = mesh.getElementNodes(1);
	EXPECT_EQ(elem1[0], 1);
	EXPECT_EQ(elem1[1], 2);
	EXPECT_EQ(elem1[2], 4);
	EXPECT_EQ(elem1[3], 5);

	const Index* elem2 = mesh.getElementNodes(2);
	EXPECT_EQ(elem2[0], 3);
	EXPECT_EQ(elem2[1], 4);
	EXPECT_EQ(elem2[2], 6);
	EXPECT_EQ(elem2[3], 7);

	const Index* elem3 = mesh.getElementNodes(3);
	EXPECT_EQ(elem3[0], 4);
	EXPECT_EQ(elem3[1], 5);
	EXPECT_EQ(elem3[2], 7);
	EXPECT_EQ(elem3[3], 8);
}

TEST(BlockMesh2D, AllElementsHaveValidNodes) {
	BlockMesh2D gen(3, 3, 0.0, 1.0, 0.0, 1.0, 1, 1);
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

TEST(BlockMesh2D, BoundaryTags_SingleElement) {
	BlockMesh2D gen(1, 1, 0.0, 1.0, 0.0, 1.0, 1, 1);
	Mesh mesh = gen.generate();

	// Single element - all faces on boundary
	// LEFT=0, RIGHT=1, BOTTOM=2, TOP=3
	EXPECT_TRUE(mesh.isOnBoundaryTag(0, 0, 0));  // LEFT
	EXPECT_TRUE(mesh.isOnBoundaryTag(0, 1, 1));  // RIGHT
	EXPECT_TRUE(mesh.isOnBoundaryTag(0, 2, 2));  // BOTTOM
	EXPECT_TRUE(mesh.isOnBoundaryTag(0, 3, 3));  // TOP
}

TEST(BlockMesh2D, BoundaryTags_2x2Mesh) {
	BlockMesh2D gen(2, 2, 0.0, 1.0, 0.0, 1.0, 1, 1);
	Mesh mesh = gen.generate();

	// Element 0 (bottom-left): LEFT and BOTTOM on boundary
	EXPECT_TRUE(mesh.isOnBoundaryTag(0, 0, 0));   // LEFT
	EXPECT_FALSE(mesh.isOnBoundary(0, 1));        // RIGHT (interior)
	EXPECT_TRUE(mesh.isOnBoundaryTag(0, 2, 2));   // BOTTOM
	EXPECT_FALSE(mesh.isOnBoundary(0, 3));        // TOP (interior)

	// Element 1 (bottom-right): RIGHT and BOTTOM on boundary
	EXPECT_FALSE(mesh.isOnBoundary(1, 0));        // LEFT (interior)
	EXPECT_TRUE(mesh.isOnBoundaryTag(1, 1, 1));   // RIGHT
	EXPECT_TRUE(mesh.isOnBoundaryTag(1, 2, 2));   // BOTTOM
	EXPECT_FALSE(mesh.isOnBoundary(1, 3));        // TOP (interior)

	// Element 2 (top-left): LEFT and TOP on boundary
	EXPECT_TRUE(mesh.isOnBoundaryTag(2, 0, 0));   // LEFT
	EXPECT_FALSE(mesh.isOnBoundary(2, 1));        // RIGHT (interior)
	EXPECT_FALSE(mesh.isOnBoundary(2, 2));        // BOTTOM (interior)
	EXPECT_TRUE(mesh.isOnBoundaryTag(2, 3, 3));   // TOP

	// Element 3 (top-right): RIGHT and TOP on boundary
	EXPECT_FALSE(mesh.isOnBoundary(3, 0));        // LEFT (interior)
	EXPECT_TRUE(mesh.isOnBoundaryTag(3, 1, 1));   // RIGHT
	EXPECT_FALSE(mesh.isOnBoundary(3, 2));        // BOTTOM (interior)
	EXPECT_TRUE(mesh.isOnBoundaryTag(3, 3, 3));   // TOP
}

TEST(BlockMesh2D, InteriorElementsNoExternalFaces) {
	// 5x5 mesh - check interior elements have no boundary faces
	BlockMesh2D gen(5, 5, 0.0, 1.0, 0.0, 1.0, 1, 1);
	Mesh mesh = gen.generate();

	// Element at (2,2) is interior
	Index interiorElem = 2 * 5 + 2;  // ele_y=2, ele_x=2

	EXPECT_FALSE(mesh.isOnBoundary(interiorElem, 0));  // LEFT
	EXPECT_FALSE(mesh.isOnBoundary(interiorElem, 1));  // RIGHT
	EXPECT_FALSE(mesh.isOnBoundary(interiorElem, 2));  // BOTTOM
	EXPECT_FALSE(mesh.isOnBoundary(interiorElem, 3));  // TOP
}

// ==================================================================
// Higher-Order Element Tests
// ==================================================================

TEST(BlockMesh2D, QuadraticElements) {
	// Quadratic elements: px=2, py=2
	BlockMesh2D gen(2, 2, 0.0, 1.0, 0.0, 1.0, 2, 2);
	Mesh mesh = gen.generate();

	EXPECT_TRUE(mesh.isValid());
	EXPECT_EQ(mesh.data.nodesPerElement, 9);
	EXPECT_EQ(mesh.data.numElements, 4);
	EXPECT_EQ(mesh.data.numNodes, 25);
	EXPECT_EQ(mesh.data.basisOrder[0], 2);
	EXPECT_EQ(mesh.data.basisOrder[1], 2);
}

TEST(BlockMesh2D, CubicElements) {
	// Cubic elements: px=3, py=3
	BlockMesh2D gen(1, 1, 0.0, 1.0, 0.0, 1.0, 3, 3);
	Mesh mesh = gen.generate();

	EXPECT_TRUE(mesh.isValid());
	EXPECT_EQ(mesh.data.nodesPerElement, 16);    // (3+1)*(3+1)
	EXPECT_EQ(mesh.data.numElements, 1);
	EXPECT_EQ(mesh.data.numNodes, 16);           // (1*3+1)*(1*3+1)
	EXPECT_EQ(mesh.data.basisOrder[0], 3);
	EXPECT_EQ(mesh.data.basisOrder[1], 3);
}

TEST(BlockMesh2D, AnisotropicOrder) {
	// Different orders in x and y: px=2, py=1
	BlockMesh2D gen(2, 2, 0.0, 1.0, 0.0, 1.0, 2, 1);
	Mesh mesh = gen.generate();

	EXPECT_TRUE(mesh.isValid());
	EXPECT_EQ(mesh.data.nodesPerElement, 6);     // (2+1)*(1+1)
	EXPECT_EQ(mesh.data.numElements, 4);         // 2*2
	EXPECT_EQ(mesh.data.numNodes, 15);           // (2*2+1)*(2*1+1)
	EXPECT_EQ(mesh.data.basisOrder[0], 2);
	EXPECT_EQ(mesh.data.basisOrder[1], 1);
}

TEST(BlockMesh2D, HigherOrderConnectivity) {
	// Verify connectivity for quadratic mesh
	BlockMesh2D gen(1, 1, 0.0, 1.0, 0.0, 1.0, 2, 2);
	Mesh mesh = gen.generate();

	// Single quadratic element should have 9 nodes
	EXPECT_EQ(mesh.data.numElements, 1);
	EXPECT_EQ(mesh.data.numNodes, 9);

	const Index* nodes = mesh.getElementNodes(0);
	// Verify all 9 nodes are unique and valid
	for (Index i = 0; i < 9; ++i) {
		EXPECT_LT(nodes[i], 9);
		for (Index j = i + 1; j < 9; ++j) {
			EXPECT_NE(nodes[i], nodes[j])
				<< "Duplicate node in element: " << nodes[i];
		}
	}
}

TEST(BlockMesh2D, HigherOrderNodeCoordinates) {
	// Single quadratic element on the unit square: the node grid must be the full
	// 3x3 fine grid (spacing 0.5), not just the 4 element-corner positions -- this is
	// a regression test for a real bug where generateNodes() ignored px/py entirely
	// and only ever placed a linear (nx+1)x(ny+1) grid, silently mislabeling every
	// mid-edge/interior node with the wrong (coarse-grid) coordinates.
	BlockMesh2D gen(1, 1, 0.0, 1.0, 0.0, 1.0, 2, 2);
	Mesh mesh = gen.generate();

	ASSERT_EQ(mesh.data.numNodes, 9);

	for (Index j = 0; j <= 2; ++j) {
		for (Index i = 0; i <= 2; ++i) {
			Index nodeID = (j * 3) + i;
			const Real* coord = mesh.getNodeCoord(nodeID);
			EXPECT_DOUBLE_EQ(coord[0], 0.5 * i) << "node " << nodeID;
			EXPECT_DOUBLE_EQ(coord[1], 0.5 * j) << "node " << nodeID;
		}
	}
}

// ==================================================================
// Mesh Validation Tests
// ==================================================================

TEST(BlockMesh2D, ValidationAfterGeneration) {
	BlockMesh2D gen(3, 4, -1.0, 2.0, 0.0, 5.0, 1, 1);
	Mesh mesh = gen.generate();

	EXPECT_TRUE(mesh.isValid());
}

TEST(BlockMesh2D, ClearMesh) {
	BlockMesh2D gen(2, 2, 0.0, 1.0, 0.0, 1.0, 1, 1);
	Mesh mesh = gen.generate();

	EXPECT_TRUE(mesh.isValid());

	mesh.clear();

	EXPECT_EQ(mesh.data.numNodes, 0);
	EXPECT_EQ(mesh.data.numElements, 0);
	EXPECT_EQ(mesh.data.xyz.size(), 0);
	EXPECT_EQ(mesh.data.ien.size(), 0);
	EXPECT_FALSE(mesh.isValid());
}
