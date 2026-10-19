#include <stdexcept>

#include <gtest/gtest.h>

#include "core/Types.hpp"
#include "io/fieldio/NodalValueSourceAdapter.hpp"
#include "mesh/generator/BlockMesh2D.hpp"

using namespace residuum;

namespace {

	// echoes the node ID back as its single component value, satisfying fem::evaluator::EvalNodalData
	struct NodeIDEchoSource {
		static constexpr Index NumComponents = 1;
		explicit NodeIDEchoSource(const mesh::Mesh&) {}
		void eval(Index nodeID, Real* out) const { out[0] = static_cast<Real>(nodeID); }
	}; // struct NodeIDEchoSource

} // namespace

class NodalValueSourceAdapterTest : public ::testing::Test {
protected:

	mesh::generator::BlockMesh2D gen{4, 4, -1.0, 1.0, -1.0, 1.0, 1, 1};
	mesh::Mesh mesh;

	void SetUp() override { mesh = gen.generate(); }

}; // class NodalValueSourceAdapterTest

TEST_F(NodalValueSourceAdapterTest, EveryMeshNodeCoordinateResolvesBackToItsOwnNodeID) {

	io::fieldio::NodalValueSourceAdapter<NodeIDEchoSource> adapter(mesh);

	for (Index nodeID = 0; nodeID < mesh.data.numNodes; ++nodeID) {
		const Real* c = mesh.getNodeCoord(nodeID);
		Real out[1];
		adapter(0.0, c, out);
		EXPECT_EQ(out[0], static_cast<Real>(nodeID));
	}

}

TEST_F(NodalValueSourceAdapterTest, ThrowsForACoordinateNotOnTheMesh) {

	io::fieldio::NodalValueSourceAdapter<NodeIDEchoSource> adapter(mesh);

	Real offMesh[2] = {123.456, -789.012};
	Real out[1];
	EXPECT_THROW(adapter(0.0, offMesh, out), std::runtime_error);

}
