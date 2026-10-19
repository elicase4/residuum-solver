#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "core/Types.hpp"
#include "fem/basis/LagrangeQuad.hpp"
#include "fem/boundary/BoundaryCategory.hpp"
#include "fem/boundary/BoundaryCondition.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/boundary/EssentialBoundaryValueLowering.hpp"
#include "mesh/generator/BlockMesh2D.hpp"
#include "topology/TopologicalDOF.hpp"

using namespace residuum;

namespace {

	// a plain position-dependent boundary function, satisfying fem::boundary::BoundaryFunction
	// without pulling in any equation-specific type
	struct LinearBoundaryFunction {
		static constexpr Index NumComponents = 1;
		Real a, b;
		void eval(Real, const Real* x, Real* value) const { value[0] = a*x[0] + b*x[1]; }
	}; // struct LinearBoundaryFunction

} // namespace

// LEFT=0 (x=x0), RIGHT=1 (x=x1), BOTTOM=2 (y=y0), TOP=3 (y=y1) -- BlockMesh2D's own convention.
class EssentialBoundaryValueLoweringTest : public ::testing::Test {
protected:

	static constexpr Index numDOFs = 1;
	static constexpr Real x0 = -1.0, x1 = 1.0, y0 = -1.0, y1 = 1.0;
	static constexpr Index nx = 4, ny = 4;
	static constexpr Real a = 1.0, b = 2.0;
	static constexpr Real sentinel = -999.0;

	fem::basis::LagrangeQuad basis{1, 1};
	mesh::generator::BlockMesh2D gen{nx, ny, x0, x1, y0, y1, 1, 1};
	mesh::Mesh mesh;
	fem::boundary::EssentialBoundaryRegistry bcRegistry;
	std::unique_ptr<topology::TopologicalDOF<numDOFs>> topoDOF;

	void SetUp() override {

		mesh = gen.generate();

		auto bc = std::make_shared<fem::boundary::BoundaryCondition<LinearBoundaryFunction>>();
		bc->tag = 0;
		bc->componentType[0] = fem::boundary::BCCategory::Essential;
		bc->function = LinearBoundaryFunction{a, b};
		bcRegistry.registerBC<LinearBoundaryFunction>(bc);

		topoDOF = std::make_unique<topology::TopologicalDOF<numDOFs>>(mesh, fem::dof::DOFOrdering::Interleaved);
		topoDOF->buildConstraints(basis, bcRegistry);

	}

}; // class EssentialBoundaryValueLoweringTest

TEST_F(EssentialBoundaryValueLoweringTest, LowersEveryConstrainedDOFToTheBoundaryFunctionsValue) {

	std::vector<Real> constrainedValues(topoDOF->numGlobalDOFs(), sentinel);
	fem::boundary::lowerEssentialBoundaryValues<numDOFs>(mesh, *topoDOF, bcRegistry, 0.0, constrainedValues.data());

	for (Index nodeID = 0; nodeID < mesh.data.numNodes; ++nodeID) {

		const Real* c = mesh.getNodeCoord(nodeID);
		Index tdof = topoDOF->getNodeDOF(nodeID, 0);

		if (c[0] == x0) {
			ASSERT_TRUE(topoDOF->isConstrained(tdof));
			EXPECT_NEAR(constrainedValues[tdof], a*c[0] + b*c[1], 1e-12);
		} else {
			EXPECT_EQ(constrainedValues[tdof], sentinel);
		}

	}

}

TEST_F(EssentialBoundaryValueLoweringTest, TimeArgumentIsForwardedToTheBoundaryFunction) {

	// LinearBoundaryFunction ignores time, so this just confirms lowerEssentialBoundaryValues
	// doesn't crash or misbehave when called again with a different time -- a real time-varying
	// boundary function would need its own test once one exists.
	std::vector<Real> constrainedValues(topoDOF->numGlobalDOFs(), sentinel);
	fem::boundary::lowerEssentialBoundaryValues<numDOFs>(mesh, *topoDOF, bcRegistry, 5.0, constrainedValues.data());

	for (Index nodeID = 0; nodeID < mesh.data.numNodes; ++nodeID) {
		const Real* c = mesh.getNodeCoord(nodeID);
		if (c[0] != x0) continue;
		Index tdof = topoDOF->getNodeDOF(nodeID, 0);
		EXPECT_NEAR(constrainedValues[tdof], a*c[0] + b*c[1], 1e-12);
	}

}
