#include <map>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "core/FEM.hpp"
#include "core/LinAlg.hpp"
#include "core/Mesh.hpp"
#include "core/Topology.hpp"
#include "core/Types.hpp"

#include "equation/heateq/HeatEquation.hpp"
#include "fem/quantity/BoundaryQuantityCombination.hpp"
#include "fem/quantity/MonitorGroup.hpp"
#include "fem/quantity/MonitorGroupImpl.hpp"
#include "fem/quantity/QuantityEvaluator.hpp"
#include "mesh/ElementFamily.hpp"
#include "mesh/generator/BlockMesh2D.hpp"

using namespace residuum;

// Exercises MonitorGroupImpl purely through the type-erased fem::quantity::MonitorGroup
// interface -- the whole point of the design (see transient-solver-design-pin memory) is that
// callers never need MonitorGroupImpl's template parameters, so these tests deliberately never
// name MonitorGroupImpl outside of construction.
class MonitorGroupImplTest : public ::testing::Test {
protected:

	// harmonic field T = a*x + b*y on a square, same setup as QuantityEvaluatorTest -- gives
	// closed-form flux integrals/averages on every boundary to check the type-erased path
	// against.
	const Real x0 = -1.0;
	const Real x1 = 1.0;
	const Real y0 = -1.0;
	const Real y1 = 1.0;
	static constexpr Index nx = 4;
	static constexpr Index ny = 4;

	static constexpr Index nsd = 2;
	static constexpr Index Px = 1;
	static constexpr Index Py = 1;
	static constexpr Index numQuadPoint = 2;

	static constexpr Real conductivity = 2.0;
	static constexpr Real a = 1.0; // dT/dx
	static constexpr Real b = 2.0; // dT/dy

	using BackendType = linalg::types::backend::CPU;
	using HeatEqBundle = equation::HeatEquation<nsd, 2, mesh::ElementFamily::Quad>;

	// LEFT=0 (x=x0), RIGHT=1 (x=x1), BOTTOM=2 (y=y0), TOP=3 (y=y1) -- BlockMesh2D's own convention.
	using Quantities = fem::quantity::QuantityForms<fem::quantity::ReducedQuantity<HeatEqBundle::HeatFluxIntegrand, fem::quantity::Reduction::Integral>>;
	using MonitorGroupT = fem::quantity::MonitorGroupImpl<BackendType, HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPBdy, HeatEqBundle::ConductivityModelBdy, Quantities, HeatEqBundle::QuadratureBoundaryType>;

	HeatEqBundle::Basis basis{Px, Py};
	HeatEqBundle::QuadratureVolumeType quadVol{numQuadPoint, numQuadPoint};
	HeatEqBundle::QuadratureBoundaryType quadBdy{numQuadPoint};
	HeatEqBundle::EvalEle evalEle{basis};

	mesh::generator::BlockMesh2D gen{nx, ny, x0, x1, y0, y1, Px, Py};
	mesh::Mesh mesh;
	std::unique_ptr<topology::TopologicalDOF<HeatEqBundle::NumDOFs>> topoDOF;

	fem::boundary::EssentialBoundaryRegistry essentialBCs;
	HeatEqBundle::ConductivityModelBdy conductivityModelBdy;
	std::unique_ptr<linalg::types::Vector<Real, BackendType>> U;

	Real fluxLeft, fluxRight, fluxBottom, fluxTop;

	void SetUp() override {

		conductivityModelBdy.setConstant(conductivity);

		mesh = gen.generate();
		topoDOF = std::make_unique<topology::TopologicalDOF<HeatEqBundle::NumDOFs>>(mesh, fem::dof::DOFOrdering::Interleaved);
		topoDOF->buildConstraints(basis, essentialBCs);

		U = std::make_unique<linalg::types::Vector<Real, BackendType>>(fem::assembly::Assembler<BackendType>::createVector<HeatEqBundle::NumDOFs>(mesh, *topoDOF));

		for (Index nodeID = 0; nodeID < mesh.data.numNodes; ++nodeID) {
			const Real* c = mesh.getNodeCoord(nodeID);
			Index tdof = topoDOF->getNodeDOF(nodeID, 0);
			Index adof = topoDOF->toAlgebraic(tdof);
			U->data()[adof] = a * c[0] + b * c[1];
		}

		fluxLeft = conductivity * a * (y1 - y0);
		fluxRight = -conductivity * a * (y1 - y0);
		fluxBottom = conductivity * b * (x1 - x0);
		fluxTop = -conductivity * b * (x1 - x0);

	}

	std::unique_ptr<fem::quantity::MonitorGroup> makeGroup() {

		Quantities quantities{typename HeatEqBundle::HeatFluxIntegrand{}};
		return std::make_unique<MonitorGroupT>(mesh, *topoDOF, essentialBCs, conductivityModelBdy, quantities, evalEle, quadBdy, *U);

	}

}; // class MonitorGroupImplTest

TEST_F(MonitorGroupImplTest, SingleTermCombinationMatchesAnalyticFluxIntegral) {

	auto group = makeGroup();

	group->registerTag(0);
	const Index combo = group->addCombination();
	group->addTerm(combo, 0, 1.0);

	bool sunk = false;
	group->evaluate(0.0, [&](Index idx, const Real* value) {
		EXPECT_EQ(idx, combo);
		EXPECT_NEAR(value[0], fluxLeft, 1e-10);
		sunk = true;
	});

	EXPECT_TRUE(sunk);

}

// Mirrors QuantityEvaluatorTest's BoundaryQuantityRegistryAndCombinationMatchAnalytic, but driven
// entirely through registerTag/addCombination/addTerm/evaluate -- i.e. through exactly the
// interface HeatProblem's monitor dispatch uses, never touching BoundaryQuantityRegistry or
// BoundaryQuantityCombination directly.
TEST_F(MonitorGroupImplTest, MultipleCombinationsShareOneMeshPassAndMatchAnalytic) {

	auto group = makeGroup();

	for (Int tag = 0; tag < 4; ++tag) group->registerTag(tag);

	const Index boundary0Alone = group->addCombination();
	group->addTerm(boundary0Alone, 0, 1.0);

	const Index sum1and2 = group->addCombination();
	group->addTerm(sum1and2, 1, 1.0);
	group->addTerm(sum1and2, 2, 1.0);

	const Index diff1minus3 = group->addCombination();
	group->addTerm(diff1minus3, 1, 1.0);
	group->addTerm(diff1minus3, 3, -1.0);

	const Index sumAll = group->addCombination();
	for (Int tag = 0; tag < 4; ++tag) group->addTerm(sumAll, tag, 1.0);

	// combination indices are handed back in creation order, starting at 0
	EXPECT_EQ(boundary0Alone, 0);
	EXPECT_EQ(sum1and2, 1);
	EXPECT_EQ(diff1minus3, 2);
	EXPECT_EQ(sumAll, 3);

	std::map<Index, Real> captured;
	group->evaluate(0.0, [&](Index idx, const Real* value) { captured[idx] = value[0]; });

	ASSERT_EQ(captured.size(), 4u);
	EXPECT_NEAR(captured[boundary0Alone], fluxLeft, 1e-10);
	EXPECT_NEAR(captured[sum1and2], fluxRight + fluxBottom, 1e-10);
	EXPECT_NEAR(captured[diff1minus3], fluxRight - fluxTop, 1e-10);
	// harmonic field, zero source -- divergence theorem says the total flux out of all four
	// boundaries is exactly zero.
	EXPECT_NEAR(captured[sumAll], 0.0, 1e-9);

}

TEST_F(MonitorGroupImplTest, EvaluateIsRepeatableAndDeterministic) {

	auto group = makeGroup();

	group->registerTag(0);
	const Index combo = group->addCombination();
	group->addTerm(combo, 0, 1.0);

	Real first = 0.0, second = 0.0;
	group->evaluate(0.0, [&](Index, const Real* value) { first = value[0]; });
	group->evaluate(1.0, [&](Index, const Real* value) { second = value[0]; }); // time is unused by a steady quantity

	EXPECT_NEAR(first, fluxLeft, 1e-10);
	EXPECT_NEAR(second, fluxLeft, 1e-10);

}
