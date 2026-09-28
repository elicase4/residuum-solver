#include <cmath>
#include <gtest/gtest.h>

#include "core/FEM.hpp"
#include "core/LinAlg.hpp"
#include "core/Mesh.hpp"
#include "core/Topology.hpp"
#include "core/Types.hpp"

#include "equation/heateq/HeatEquation.hpp"
#include "fem/quantity/BoundaryQuantityCombination.hpp"
#include "fem/quantity/QuantityEvaluator.hpp"
#include "mesh/ElementFamily.hpp"
#include "mesh/generator/BlockMesh2D.hpp"

using namespace residuum;

namespace {

	struct BoundaryLengthQuantity {
		static constexpr Index NumComponents = 1;
		void computeElementLevelValue(const auto& qp, const Real*, Real* out) const {
			using QP = std::decay_t<decltype(qp)>;
			Real normalMag = 0;
			for (Index sD = 0; sD < QP::SpatialDim; ++sD) normalMag += qp.normal[sD] * qp.normal[sD];
			out[0] = std::sqrt(normalMag) * qp.w;
		}
	};

} // namespace

class QuantityEvaluatorTest : public ::testing::Test {
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

	HeatEqBundle::Basis basis{Px, Py};
	HeatEqBundle::QuadratureVolumeType quadVol{numQuadPoint, numQuadPoint};
	HeatEqBundle::QuadratureBoundaryType quadBdy{numQuadPoint};
	HeatEqBundle::EvalEle evalEle{basis};

	mesh::generator::BlockMesh2D gen{nx, ny, x0, x1, y0, y1, Px, Py};
	mesh::Mesh mesh;
	std::unique_ptr<topology::TopologicalDOF<HeatEqBundle::NumDOFs>> topoDOF;

	// no essential BCs registered
	fem::boundary::EssentialBoundaryRegistry essentialBCs;

	HeatEqBundle::ConductivityModelBdy conductivityModelBdy;

	std::unique_ptr<linalg::types::Vector<Real, BackendType>> U;

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

	}

};

// LEFT boundary (tag 0, x = x0): outward normal (-1, 0), length (y1 - y0). flux = -K*gradT.n = -k*(a,b).(-1,0) = k*a, uniform along the whole boundary.
// evaluateBoundary was removed (fully subsumed by evaluateBoundaryRegistry with one registered
// tag) -- single-tag cases below go through the registry API now.
TEST_F(QuantityEvaluatorTest, HeatFluxAverageOnLeftBoundaryMatchesAnalytic) {

	static constexpr Int boundaryTag = 0;
	const Real expectedAverage = conductivity * a;

	using Quantities = fem::quantity::QuantityForms<fem::quantity::ReducedQuantity<HeatEqBundle::HeatFluxIntegrand, fem::quantity::Reduction::Average>>;
	static_assert(Quantities::TotalComponents <= equation::heateq::quantity::kMaxQuantityComponents);

	Quantities quantities{typename HeatEqBundle::HeatFluxIntegrand{}};

	fem::quantity::BoundaryQuantityRegistry<Quantities> registry;
	registry.registerTag(boundaryTag);

	fem::quantity::QuantityEvaluator<BackendType>::evaluateBoundaryRegistry<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPBdy, HeatEqBundle::ConductivityModelBdy, Quantities, HeatEqBundle::QuadratureBoundaryType>(mesh, *topoDOF, essentialBCs, 0.0, conductivityModelBdy, quantities, evalEle, quadBdy, *U, registry);

	EXPECT_NEAR(registry.result(boundaryTag)[Quantities::offset<0>()], expectedAverage, 1e-10);

}

// Two independent quantities, different types, same target, in ONE mesh pass -- the flux integral plus a trivial "boundary length" quantity that should come out to exactly (y1 - y0).
TEST_F(QuantityEvaluatorTest, MultipleQuantitiesInOnePassAreIndependentlyCorrect) {

	static constexpr Int boundaryTag = 0;
	const Real expectedFluxIntegral = conductivity * a * (y1 - y0);
	const Real expectedLength = y1 - y0;

	using Quantities = fem::quantity::QuantityForms<
		fem::quantity::ReducedQuantity<HeatEqBundle::HeatFluxIntegrand, fem::quantity::Reduction::Integral>,
		fem::quantity::ReducedQuantity<BoundaryLengthQuantity, fem::quantity::Reduction::Integral>
	>;
	static_assert(Quantities::TotalComponents <= equation::heateq::quantity::kMaxQuantityComponents);

	Quantities quantities{typename HeatEqBundle::HeatFluxIntegrand{}, BoundaryLengthQuantity{}};

	fem::quantity::BoundaryQuantityRegistry<Quantities> registry;
	registry.registerTag(boundaryTag);

	fem::quantity::QuantityEvaluator<BackendType>::evaluateBoundaryRegistry<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPBdy, HeatEqBundle::ConductivityModelBdy, Quantities, HeatEqBundle::QuadratureBoundaryType>(mesh, *topoDOF, essentialBCs, 0.0, conductivityModelBdy, quantities, evalEle, quadBdy, *U, registry);

	EXPECT_NEAR(registry.result(boundaryTag)[Quantities::offset<0>()], expectedFluxIntegral, 1e-10);
	EXPECT_NEAR(registry.result(boundaryTag)[Quantities::offset<1>()], expectedLength, 1e-10);

}

// Registers a real Dirichlet BC (matching T exactly) on every boundary, including the flux
// target itself -- corner nodes shared between LEFT and its neighbors are now CONSTRAINED, so
// this only reproduces the exact analytic flux if the constrained-DOF gather is actually correct.
TEST_F(QuantityEvaluatorTest, HeatFluxIntegralIsCorrectWhenTargetBoundaryIsAlsoDirichletConstrained) {

	static constexpr Int boundaryTag = 0;
	const Real expectedIntegral = conductivity * a * (y1 - y0);

	auto g = [](Real, const Real* x, Real* out) { out[0] = a * x[0] + b * x[1]; };
	using DirichletT = HeatEqBundle::DirichletExpression<decltype(g)>;

	fem::boundary::EssentialBoundaryRegistry constrainedBCs;
	for (Int boundaryID = 0; boundaryID < 4; ++boundaryID) {
		auto bc = std::shared_ptr<fem::boundary::BoundaryCondition<DirichletT>>(new fem::boundary::BoundaryCondition<DirichletT>{boundaryID, {fem::boundary::BCCategory::Essential}, DirichletT{g}});
		constrainedBCs.registerBC<DirichletT>(bc);
	}

	auto constrainedTopoDOF = std::make_unique<topology::TopologicalDOF<HeatEqBundle::NumDOFs>>(mesh, fem::dof::DOFOrdering::Interleaved);
	constrainedTopoDOF->buildConstraints(basis, constrainedBCs);

	auto constrainedU = fem::assembly::Assembler<BackendType>::createVector<HeatEqBundle::NumDOFs>(mesh, *constrainedTopoDOF);
	for (Index nodeID = 0; nodeID < mesh.data.numNodes; ++nodeID) {
		Index tdof = constrainedTopoDOF->getNodeDOF(nodeID, 0);
		if (constrainedTopoDOF->isConstrained(tdof)) continue; // left for constrainedBCs to supply
		const Real* c = mesh.getNodeCoord(nodeID);
		Index adof = constrainedTopoDOF->toAlgebraic(tdof);
		constrainedU.data()[adof] = a * c[0] + b * c[1];
	}

	using Quantities = fem::quantity::QuantityForms<fem::quantity::ReducedQuantity<HeatEqBundle::HeatFluxIntegrand, fem::quantity::Reduction::Integral>>;
	Quantities quantities{typename HeatEqBundle::HeatFluxIntegrand{}};

	fem::quantity::BoundaryQuantityRegistry<Quantities> registry;
	registry.registerTag(boundaryTag);

	fem::quantity::QuantityEvaluator<BackendType>::evaluateBoundaryRegistry<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPBdy, HeatEqBundle::ConductivityModelBdy, Quantities, HeatEqBundle::QuadratureBoundaryType>(mesh, *constrainedTopoDOF, constrainedBCs, 0.0, conductivityModelBdy, quantities, evalEle, quadBdy, constrainedU, registry);

	EXPECT_NEAR(registry.result(boundaryTag)[Quantities::offset<0>()], expectedIntegral, 1e-10);

}

// LEFT=0 (x=x0), RIGHT=1 (x=x1), BOTTOM=2 (y=y0), TOP=3 (y=y1) -- BlockMesh2D's own convention.
// One registry pass computes all four tags' flux integrals; combinations are checked against
// their closed-form values, including the physically meaningful one: since T is harmonic (zero
// source), the total flux out of all four boundaries must be exactly zero (divergence theorem).
TEST_F(QuantityEvaluatorTest, BoundaryQuantityRegistryAndCombinationMatchAnalytic) {

	using Quantities = fem::quantity::QuantityForms<fem::quantity::ReducedQuantity<HeatEqBundle::HeatFluxIntegrand, fem::quantity::Reduction::Integral>>;
	Quantities quantities{typename HeatEqBundle::HeatFluxIntegrand{}};

	fem::quantity::BoundaryQuantityRegistry<Quantities> registry;
	for (Int tag = 0; tag < 4; ++tag) registry.registerTag(tag);

	fem::quantity::QuantityEvaluator<BackendType>::evaluateBoundaryRegistry<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPBdy, HeatEqBundle::ConductivityModelBdy, Quantities, HeatEqBundle::QuadratureBoundaryType>(mesh, *topoDOF, essentialBCs, 0.0, conductivityModelBdy, quantities, evalEle, quadBdy, *U, registry);

	const Real fluxLeft = conductivity * a * (y1 - y0);
	const Real fluxRight = -conductivity * a * (y1 - y0);
	const Real fluxBottom = conductivity * b * (x1 - x0);
	const Real fluxTop = -conductivity * b * (x1 - x0);

	EXPECT_NEAR(registry.result(0)[0], fluxLeft, 1e-10);
	EXPECT_NEAR(registry.result(1)[0], fluxRight, 1e-10);
	EXPECT_NEAR(registry.result(2)[0], fluxBottom, 1e-10);
	EXPECT_NEAR(registry.result(3)[0], fluxTop, 1e-10);

	Real out[HeatEqBundle::HeatFluxIntegrand::NumComponents];

	fem::quantity::BoundaryQuantityCombination<Quantities> boundary0Alone;
	boundary0Alone.addTerm(0, 1.0);
	boundary0Alone.evaluate(registry, out);
	EXPECT_NEAR(out[0], fluxLeft, 1e-10);

	fem::quantity::BoundaryQuantityCombination<Quantities> sum1and2;
	sum1and2.addTerm(1, 1.0);
	sum1and2.addTerm(2, 1.0);
	sum1and2.evaluate(registry, out);
	EXPECT_NEAR(out[0], fluxRight + fluxBottom, 1e-10);

	fem::quantity::BoundaryQuantityCombination<Quantities> diff1minus3;
	diff1minus3.addTerm(1, 1.0);
	diff1minus3.addTerm(3, -1.0);
	diff1minus3.evaluate(registry, out);
	EXPECT_NEAR(out[0], fluxRight - fluxTop, 1e-10);

	fem::quantity::BoundaryQuantityCombination<Quantities> sumAll;
	for (Int tag = 0; tag < 4; ++tag) sumAll.addTerm(tag, 1.0);
	sumAll.evaluate(registry, out);
	EXPECT_NEAR(out[0], 0.0, 1e-9);

}

TEST_F(QuantityEvaluatorTest, BoundaryQuantityRegistryThrowsOnUnregisteredTag) {

	using Quantities = fem::quantity::QuantityForms<fem::quantity::ReducedQuantity<HeatEqBundle::HeatFluxIntegrand, fem::quantity::Reduction::Integral>>;
	fem::quantity::BoundaryQuantityRegistry<Quantities> registry;
	registry.registerTag(0);

	EXPECT_THROW(registry.result(1), std::runtime_error);

}
