#include <memory>
#include <utility>

#include <gtest/gtest.h>

#include "core/Types.hpp"

#include "application/heateq/config/HeatConfig.hpp"
#include "application/heateq/problem/HeatProblem.hpp"

#include "equation/heateq/HeatEquation.hpp"

#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/dof/DOFOrdering.hpp"

#include "linalg/types/backend/CPU.hpp"

#include "mesh/ElementFamily.hpp"
#include "mesh/generator/BlockMesh2D.hpp"

#include "solver/config/DriverConfig.hpp"
#include "solver/stage/BackwardEulerStage.hpp"
#include "solver/stage/TransientCapableStage.hpp"

#include "topology/TopologicalDOF.hpp"

using namespace residuum;

namespace {

	using BackendType = linalg::types::backend::CPU;
	using HeatEqBundle = equation::HeatEquation<2, 2, mesh::ElementFamily::Quad>;
	using HeatProblemT = application::heateq::problem::HeatProblem<BackendType, HeatEqBundle>;
	using BEStageT = solver::stage::BackwardEulerStage<HeatProblemT>;

} // namespace

static_assert(solver::stage::Stage<BEStageT>);
static_assert(solver::stage::TransientCapableStage<BEStageT>);

class BackwardEulerStageTest : public ::testing::Test {
protected:

	const Real x0 = -1.0, x1 = 1.0, y0 = -1.0, y1 = 1.0;
	static constexpr Index nx = 4, ny = 4;
	static constexpr Index Px = 1, Py = 1;

	static constexpr Real conductivity = 2.0;
	static constexpr Real density = 1.0;
	static constexpr Real specificHeat = 1.0;
	static constexpr Real a = 1.0; // dT/dx
	static constexpr Real b = 2.0; // dT/dy

	application::heateq::config::HeatConfig makeConfig() const {

		namespace hconfig = application::heateq::config;
		namespace sconfig = solver::config;

		hconfig::HeatConfig cfg;

		cfg.discretization.quadrature.xi = 2;
		cfg.discretization.quadrature.eta = 2;
		cfg.discretization.dofOrdering = fem::dof::DOFOrdering::Interleaved;

		sconfig::LinearSolverConfig linSolver;
		linSolver.type = sconfig::LinearSolverConfig::Type::CG;
		linSolver.operatorType = sconfig::LinearSolverConfig::OperatorType::CSR;
		linSolver.tolerance = 1e-12;
		linSolver.maxIterations = 2000;
		cfg.solver.linear = linSolver;

		cfg.solver.driver.type = sconfig::DriverConfig::Type::Transient;

		sconfig::TimeStepperConfig ts;
		ts.type = sconfig::TimeStepperConfig::Type::BackwardEuler;
		ts.t0 = 0.0;
		ts.tf = 1.0;
		ts.stepSize.mode = sconfig::TimeStepSizeConfig::Mode::Constant;
		ts.stepSize.dt = 0.1;
		ts.linearSolver = linSolver;
		cfg.solver.timestepper = ts;

		cfg.density = hconfig::ScalarMaterialPropertyConfig{density, "kg/m^3"};
		cfg.specificHeat = hconfig::SpecificHeatConfig{hconfig::SpecificHeatConfig::Type::Constant, specificHeat, "", "", "J/(kg*K)"};

		cfg.conductivity.type = hconfig::ConductivityConfig::Type::Constant;
		cfg.conductivity.value = conductivity;
		cfg.conductivity.unit = "W/(m*K)";

		cfg.source.type = hconfig::SourceConfig::Type::VolumetricHeatSource;
		cfg.source.read.mode = sconfig::NodalFieldReadConfig::Mode::Expression;
		cfg.source.read.expression = "0.0";
		cfg.source.read.unit = "W/m^3";

		// steady-state IC, see class comment
		cfg.initialCondition.read.mode = sconfig::NodalFieldReadConfig::Mode::Expression;
		cfg.initialCondition.read.expression = "1.0*x + 2.0*y";
		cfg.initialCondition.read.unit = "K";

		for (Int tag = 0; tag < 4; ++tag) {
			hconfig::BoundaryConditionConfig bc;
			bc.boundaryID = tag;
			bc.type = hconfig::BoundaryConditionConfig::Type::Value;
			bc.mode = sconfig::NodalFieldReadConfig::Mode::Expression;
			bc.expression = "1.0*x + 2.0*y";
			bc.unit = "K";
			bc.forms = {hconfig::BoundaryConditionConfig::Form::ValueBC};
			bc.model = hconfig::ConductivityConfig::Type::Constant;
			cfg.boundaryConditions.push_back(bc);
		}

		return cfg;

	}

	std::unique_ptr<HeatProblemT> makeProblem() const {

		const auto cfg = makeConfig();

		HeatEqBundle::Basis basis{Px, Py};
		HeatEqBundle::QuadratureVolumeType quadVol{2, 2};
		HeatEqBundle::QuadratureBoundaryType quadBdy{2};

		mesh::generator::BlockMesh2D gen{nx, ny, x0, x1, y0, y1, Px, Py};
		mesh::Mesh mesh = gen.generate();

		return std::make_unique<HeatProblemT>(cfg, std::move(mesh), std::move(basis), std::move(quadVol), std::move(quadBdy));

	}

	struct VerificationTopology {

		mesh::Mesh mesh;
		fem::boundary::EssentialBoundaryRegistry bcs;
		std::unique_ptr<topology::TopologicalDOF<HeatEqBundle::NumDOFs>> topoDOF;

		VerificationTopology(Index nx, Index ny, Real x0, Real x1, Real y0, Real y1, Index Px, Index Py, Real a, Real b) {

			mesh::generator::BlockMesh2D gen{nx, ny, x0, x1, y0, y1, Px, Py};
			mesh = gen.generate();

			auto g = [a, b](Real, const Real* xyz, Real* out) { out[0] = a * xyz[0] + b * xyz[1]; };
			using DirichletT = HeatEqBundle::DirichletExpression<decltype(g)>;

			for (Int tag = 0; tag < 4; ++tag) {
				auto bc = std::shared_ptr<fem::boundary::BoundaryCondition<DirichletT>>(new fem::boundary::BoundaryCondition<DirichletT>{tag, {fem::boundary::BCCategory::Essential}, DirichletT{g}});
				bcs.registerBC<DirichletT>(bc);
			}

			topoDOF = std::make_unique<topology::TopologicalDOF<HeatEqBundle::NumDOFs>>(mesh, fem::dof::DOFOrdering::Interleaved);

			HeatEqBundle::Basis basis{Px, Py};
			topoDOF->buildConstraints(basis, bcs);

		}

	}; // struct VerificationTopology

	// Regression config for the TangentMassForm dt-scaling bug: BackwardEulerStage::setDt()
	// used to wrap TangentMassForms in an extra ScaledForm(1/dt) on top of the 1/dt already
	// baked into Udot, corrupting the mass tangent by a factor of dt whenever dt != 1 -- so this
	// needs (a) temperature-dependent specific heat, so the tangent term is actually nonzero,
	// and (b) dt far from 1, so a wrong scale factor is far from a no-op.
	application::heateq::config::HeatConfig makeNonlinearConfig(Real dt) const {

		namespace hconfig = application::heateq::config;
		namespace sconfig = solver::config;

		hconfig::HeatConfig cfg;

		cfg.discretization.quadrature.xi = 2;
		cfg.discretization.quadrature.eta = 2;
		cfg.discretization.dofOrdering = fem::dof::DOFOrdering::Interleaved;

		sconfig::LinearSolverConfig linSolver;
		linSolver.type = sconfig::LinearSolverConfig::Type::CG;
		linSolver.operatorType = sconfig::LinearSolverConfig::OperatorType::CSR;
		linSolver.tolerance = 1e-13;
		linSolver.maxIterations = 2000;

		sconfig::NonlinearSolverConfig nl;
		nl.type = sconfig::NonlinearSolverConfig::Type::Newton;
		nl.absoluteTolerance = 1e-14;
		nl.relativeTolerance = 1e-9;
		nl.maxIterations = 10; // quadratic Newton needs a handful; the dt-scaling bug needs ~40+
		nl.linearSolver = linSolver;
		cfg.solver.nonlinear = nl;

		cfg.solver.driver.type = sconfig::DriverConfig::Type::Transient;

		sconfig::TimeStepperConfig ts;
		ts.type = sconfig::TimeStepperConfig::Type::BackwardEuler;
		ts.t0 = 0.0;
		ts.tf = dt;
		ts.stepSize.mode = sconfig::TimeStepSizeConfig::Mode::Constant;
		ts.stepSize.dt = dt;
		ts.nonlinearSolver = nl;
		ts.linearSolver = linSolver;
		cfg.solver.timestepper = ts;

		cfg.density = hconfig::ScalarMaterialPropertyConfig{density, "kg/m^3"};
		cfg.specificHeat = hconfig::SpecificHeatConfig{hconfig::SpecificHeatConfig::Type::TemperatureDependent, 0.0, "1.0 + 0.1*T", "0.1", "J/(kg*K)"};

		cfg.conductivity.type = hconfig::ConductivityConfig::Type::Constant;
		cfg.conductivity.value = conductivity;
		cfg.conductivity.unit = "W/(m*K)";

		cfg.source.type = hconfig::SourceConfig::Type::VolumetricHeatSource;
		cfg.source.read.mode = sconfig::NodalFieldReadConfig::Mode::Expression;
		cfg.source.read.expression = "0.0";
		cfg.source.read.unit = "W/m^3";

		// starts uniformly cold; the boundary jump to a non-uniform field forces several
		// genuine Newton iterations, unlike the steady-state IC/BC trick used elsewhere in
		// this file (matched IC/BC there gives ~zero residual from iteration 0, which would
		// never touch the Jacobian's mass tangent at all)
		cfg.initialCondition.read.mode = sconfig::NodalFieldReadConfig::Mode::Expression;
		cfg.initialCondition.read.expression = "0.0";
		cfg.initialCondition.read.unit = "K";

		for (Int tag = 0; tag < 4; ++tag) {
			hconfig::BoundaryConditionConfig bc;
			bc.boundaryID = tag;
			bc.type = hconfig::BoundaryConditionConfig::Type::Value;
			bc.mode = sconfig::NodalFieldReadConfig::Mode::Expression;
			bc.expression = "1.0*x + 2.0*y";
			bc.unit = "K";
			bc.forms = {hconfig::BoundaryConditionConfig::Form::ValueBC};
			bc.model = hconfig::ConductivityConfig::Type::Constant;
			cfg.boundaryConditions.push_back(bc);
		}

		return cfg;

	}

	std::unique_ptr<HeatProblemT> makeNonlinearProblem(Real dt) const {

		const auto cfg = makeNonlinearConfig(dt);

		HeatEqBundle::Basis basis{Px, Py};
		HeatEqBundle::QuadratureVolumeType quadVol{2, 2};
		HeatEqBundle::QuadratureBoundaryType quadBdy{2};

		mesh::generator::BlockMesh2D gen{nx, ny, x0, x1, y0, y1, Px, Py};
		mesh::Mesh mesh = gen.generate();

		return std::make_unique<HeatProblemT>(cfg, std::move(mesh), std::move(basis), std::move(quadVol), std::move(quadBdy));

	}

	void expectSolutionMatchesAnalyticField(const HeatProblemT& problem) const {

		VerificationTopology verify(nx, ny, x0, x1, y0, y1, Px, Py, a, b);
		const Real* solutionData = problem.solution().data();

		Index numFreeChecked = 0;
		for (Index nodeID = 0; nodeID < verify.mesh.data.numNodes; ++nodeID) {

			const Index tdof = verify.topoDOF->getNodeDOF(nodeID, 0);
			if (verify.topoDOF->isConstrained(tdof)) continue;

			const Index adof = verify.topoDOF->toAlgebraic(tdof);
			const Real* c = verify.mesh.getNodeCoord(nodeID);
			const Real expected = a * c[0] + b * c[1];

			EXPECT_NEAR(solutionData[adof], expected, 1e-8);
			++numFreeChecked;

		}

		EXPECT_GT(numFreeChecked, 0);

	}

}; // class BackwardEulerStageTest

TEST_F(BackwardEulerStageTest, SteadyStateInitialConditionStaysExactAcrossMultipleTimesteps) {

	auto problem = makeProblem();
	BEStageT stage(*problem);

	const Real dt = 0.1;
	Real time = 0.0;

	for (Index step = 0; step < 3; ++step) {

		time += dt;

		stage.setDt(dt);
		stage.setTime(time);
		stage.assemble();

		ASSERT_TRUE(stage.solve());

		stage.advance();
		stage.onStepComplete(step + 1, time); // no-op here: no output/monitors configured

		expectSolutionMatchesAnalyticField(*problem);

	}

}

TEST_F(BackwardEulerStageTest, AdvanceCopiesCurrentSolutionIntoUPrev) {

	auto problem = makeProblem();
	BEStageT stage(*problem);

	stage.setDt(0.1);
	stage.setTime(0.1);
	stage.assemble();
	ASSERT_TRUE(stage.solve());
	stage.advance();

	ASSERT_EQ(problem->U_prev().size(), problem->U().size());
	for (Index i = 0; i < problem->U().size(); ++i) {
		EXPECT_DOUBLE_EQ(problem->U_prev().data()[i], problem->U().data()[i]);
	}

}

TEST_F(BackwardEulerStageTest, NewtonConvergesWithinFewIterationsForTemperatureDependentCapacityAtNonUnitDt) {

	const Real dt = 0.05; // far from 1, so a wrong 1/dt vs 1/dt^2 mass-tangent scale is not a no-op

	auto problem = makeNonlinearProblem(dt);
	BEStageT stage(*problem);

	stage.setDt(dt);
	stage.setTime(dt);
	stage.assemble();

	// with the correct Jacobian this is a handful of quadratic-ish Newton iterations; with the
	// TangentMassForm dt-scaling bug (see BackwardEulerStage::setDt()) convergence degrades to a
	// constant ~0.6 ratio per iteration, needing ~40+ iterations to reach the same tolerance --
	// far past this stage's configured max of 10, so this fails loudly if the bug regresses
	EXPECT_TRUE(stage.solve());

}

TEST_F(BackwardEulerStageTest, SolutionAccessorAliasesProblemU) {

	auto problem = makeProblem();
	BEStageT stage(*problem);

	stage.setDt(0.1);
	stage.setTime(0.1);
	stage.assemble();
	ASSERT_TRUE(stage.solve());

	EXPECT_EQ(&stage.solution(), &problem->U());

}
