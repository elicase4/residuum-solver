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
#include "equation/heateq/evaluator/HeatStateVector.hpp"
#include "fem/quantity/BoundaryReactionMonitorGroup.hpp"
#include "fem/quantity/MonitorGroup.hpp"
#include "fem/quantity/QuantityEvaluator.hpp"
#include "mesh/ElementFamily.hpp"
#include "mesh/generator/BlockMesh2D.hpp"

using namespace residuum;

// Exercises BoundaryReactionMonitorGroup purely through the type-erased fem::quantity::MonitorGroup
// interface, same spirit as BoundaryQuadratureMonitorGroupTest. Uses the identical harmonic field
// (T = a*x + b*y, constant conductivity, no source) as that test, letting this test reuse its
// fluxLeft/fluxRight/fluxBottom/fluxTop analytic values -- but a boundary REACTION (K*T - F at a
// boundary row, the FEM reaction-force convention, generic across PDEs) is not itself the physical
// outward flux: the weak form's boundary term is v*(k*gradT.n), and physical flux uses
// q = -k*gradT, so the reaction is exactly -1 times the outward flux computed by the strong path.
// heateq::problem::HeatProblem is the layer that knows this and negates accordingly when exposing
// a weak-method heat_flux monitor to match the strong method's sign; this test checks the generic
// evaluator honestly, against -flux*.
class BoundaryReactionMonitorGroupTest : public ::testing::Test {
protected:

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

	using VectorT = linalg::types::Vector<Real, BackendType>;
	using StateVectorT = equation::heateq::evaluator::HeatStateVector<HeatEqBundle::NumDOFs, HeatEqBundle::SpatialDim, HeatEqBundle::EvalEle::ParametricDim, VectorT>;

	static Real* selectValue(StateVectorT& s) { return s.valueBuffer(); }
	using SelectorT = decltype(&selectValue);

	// steady, no source: FormsT is stiffness alone, operand is the current solution T
	using FormsT = fem::form::FormRegistry<HeatEqBundle::DiffusionForm>;
	using PairT = fem::assembly::FormsWithOperand<FormsT, SelectorT>;
	using MonitorGroupT = fem::quantity::BoundaryReactionMonitorGroup<BackendType, HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::EvalQPBdy, HeatEqBundle::ConductivityModel, HeatEqBundle::QuadratureVolumeType, HeatEqBundle::QuadratureBoundaryType, StateVectorT, PairT>;

	HeatEqBundle::Basis basis{Px, Py};
	HeatEqBundle::QuadratureVolumeType quadVol{numQuadPoint, numQuadPoint};
	HeatEqBundle::QuadratureBoundaryType quadBdy{numQuadPoint};
	HeatEqBundle::EvalEle evalEle{basis};

	mesh::generator::BlockMesh2D gen{nx, ny, x0, x1, y0, y1, Px, Py};
	mesh::Mesh mesh;
	std::unique_ptr<topology::TopologicalDOF<HeatEqBundle::NumDOFs>> topoDOF;

	fem::boundary::EssentialBoundaryRegistry essentialBCs;
	HeatEqBundle::ConductivityModel conductivityModel;
	std::unique_ptr<VectorT> U;

	FormsT forms;
	StateVectorT state;

	Real fluxLeft, fluxRight, fluxBottom, fluxTop;

	void SetUp() override {

		conductivityModel.setConstant(conductivity);

		mesh = gen.generate();
		topoDOF = std::make_unique<topology::TopologicalDOF<HeatEqBundle::NumDOFs>>(mesh, fem::dof::DOFOrdering::Interleaved);
		topoDOF->buildConstraints(basis, essentialBCs);

		U = std::make_unique<VectorT>(fem::assembly::Assembler<BackendType>::createVector<HeatEqBundle::NumDOFs>(mesh, *topoDOF));

		for (Index nodeID = 0; nodeID < mesh.data.numNodes; ++nodeID) {
			const Real* c = mesh.getNodeCoord(nodeID);
			Index tdof = topoDOF->getNodeDOF(nodeID, 0);
			Index adof = topoDOF->toAlgebraic(tdof);
			U->data()[adof] = a * c[0] + b * c[1];
		}

		state.valueSource = U.get();

		fluxLeft = conductivity * a * (y1 - y0);
		fluxRight = -conductivity * a * (y1 - y0);
		fluxBottom = conductivity * b * (x1 - x0);
		fluxTop = -conductivity * b * (x1 - x0);

	}

	std::unique_ptr<fem::quantity::MonitorGroup> makeGroup(fem::quantity::Reduction reduction) {
		return std::make_unique<MonitorGroupT>(mesh, *topoDOF, essentialBCs, conductivityModel, evalEle, quadVol, quadBdy, state, reduction, PairT{forms, &selectValue});
	}

}; // class BoundaryReactionMonitorGroupTest

TEST_F(BoundaryReactionMonitorGroupTest, SingleTermCombinationMatchesAnalyticFluxIntegral) {

	auto group = makeGroup(fem::quantity::Reduction::Integral);

	group->registerTag(0);
	const Index combo = group->addCombination();
	group->addTerm(combo, 0, 1.0);

	bool sunk = false;
	group->evaluate(0.0, [&](Index idx, const Real* value) {
		EXPECT_EQ(idx, combo);
		EXPECT_NEAR(value[0], -fluxLeft, 1e-10);
		sunk = true;
	});

	EXPECT_TRUE(sunk);

}

TEST_F(BoundaryReactionMonitorGroupTest, MultipleCombinationsShareOneMeshPassAndMatchAnalytic) {

	auto group = makeGroup(fem::quantity::Reduction::Integral);

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

	std::map<Index, Real> captured;
	group->evaluate(0.0, [&](Index idx, const Real* value) { captured[idx] = value[0]; });

	ASSERT_EQ(captured.size(), 4u);
	EXPECT_NEAR(captured[boundary0Alone], -fluxLeft, 1e-10);
	EXPECT_NEAR(captured[sum1and2], -(fluxRight + fluxBottom), 1e-10);
	EXPECT_NEAR(captured[diff1minus3], -(fluxRight - fluxTop), 1e-10);
	// harmonic field, zero source -- divergence theorem says the total flux out of all four
	// boundaries is exactly zero, same as the strong path; sign-invariant, so this one doesn't
	// distinguish reaction from flux, but it does catch the corner-node dedup bug across tags.
	EXPECT_NEAR(captured[sumAll], 0.0, 1e-9);

}

TEST_F(BoundaryReactionMonitorGroupTest, AverageModeDividesByBoundaryMeasure) {

	auto group = makeGroup(fem::quantity::Reduction::Average);

	group->registerTag(0);
	const Index combo = group->addCombination();
	group->addTerm(combo, 0, 1.0);

	group->evaluate(0.0, [&](Index, const Real* value) {
		EXPECT_NEAR(value[0], -fluxLeft / (y1 - y0), 1e-10);
	});

}

TEST_F(BoundaryReactionMonitorGroupTest, EvaluateIsRepeatableAndDeterministic) {

	auto group = makeGroup(fem::quantity::Reduction::Integral);

	group->registerTag(0);
	const Index combo = group->addCombination();
	group->addTerm(combo, 0, 1.0);

	Real first = 0.0, second = 0.0;
	group->evaluate(0.0, [&](Index, const Real* value) { first = value[0]; });
	group->evaluate(1.0, [&](Index, const Real* value) { second = value[0]; }); // time is unused by a steady quantity

	EXPECT_NEAR(first, -fluxLeft, 1e-10);
	EXPECT_NEAR(second, -fluxLeft, 1e-10);

}
