#include <cmath>
#include <limits>
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
#include "solver/config/LoggingConfig.hpp"
#include "solver/driver/Transient.hpp"
#include "solver/stage/BackwardEulerStage.hpp"
#include "solver/timestepper/TimeStepperFactory.hpp"

#include "topology/TopologicalDOF.hpp"

using namespace residuum;

// End-to-end transient solve through the REAL production path -- HeatProblem +
// BackwardEulerStage + TimeStepperFactory::makeTimeStepperRunner + solver::driver::Transient --
// as opposed to UnitTest_Solver_Stage_BackwardEulerStage, which drives Stage methods by hand and
// never touches the timestepper/driver layer. This is what HeatDispatcher itself calls.
//
// Physical scenario: 1D transient conduction, folded into a 2D mesh by insulating (leaving
// unconstrained -- natural/zero-flux by default) the top and bottom edges so the field stays
// uniform in y. Domain x in [0, L], Dirichlet T=0 at both ends, zero source, IC T(x,0) =
// sin(pi*x/L) -- a pure eigenmode of the Laplacian, so the exact solution is the classic
// separation-of-variables decay:
//
//   T(x,t) = sin(pi*x/L) * exp(-alpha*(pi/L)^2*t),  alpha = k/(rho*cp)
//
// With k = rho = cp = 1 and L = 1, alpha = 1 and the decay rate is exp(-pi^2*t).
class HeatEquationTransientBackwardEulerTest : public ::testing::Test {
protected:

	const Real x0 = 0.0, x1 = 1.0, y0 = 0.0, y1 = 1.0;
	static constexpr Index nx = 20, ny = 2;
	static constexpr Index Px = 1, Py = 1;

	static constexpr Real conductivity = 1.0;
	static constexpr Real density = 1.0;
	static constexpr Real specificHeat = 1.0;
	static constexpr Real L = 1.0;
	static constexpr Real dt = 0.01;
	static constexpr Real tf = 0.1;

	using BackendType = linalg::types::backend::CPU;
	using HeatEqBundle = equation::HeatEquation<2, 2, mesh::ElementFamily::Quad>;
	using HeatProblemT = application::heateq::problem::HeatProblem<BackendType, HeatEqBundle>;
	using BEStageT = solver::stage::BackwardEulerStage<HeatProblemT>;

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
		ts.tf = tf;
		ts.stepSize.mode = sconfig::TimeStepSizeConfig::Mode::Constant;
		ts.stepSize.dt = dt;
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

		cfg.initialCondition.read.mode = sconfig::NodalFieldReadConfig::Mode::Expression;
		cfg.initialCondition.read.expression = "sin(pi*x)";
		cfg.initialCondition.read.unit = "K";

		// only LEFT (0) and RIGHT (1) are constrained -- BOTTOM (2) and TOP (3) are left
		// unregistered, which is exactly the natural/zero-flux (insulated) condition, keeping
		// the field uniform in y so this reduces to the 1D analytic case above.
		for (Int tag : {0, 1}) {
			hconfig::BoundaryConditionConfig bc;
			bc.boundaryID = tag;
			bc.type = hconfig::BoundaryConditionConfig::Type::Value;
			bc.mode = sconfig::NodalFieldReadConfig::Mode::Expression;
			bc.expression = "0.0";
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

	// finds the free-DOF algebraic index closest to (x, y) in an independently rebuilt,
	// identically-constrained topology -- see the same pattern in
	// UnitTest_Application_Heateq_Problem_HeatProblem / UnitTest_Solver_Stage_*Stage.
	struct VerificationTopology {

		mesh::Mesh mesh;
		fem::boundary::EssentialBoundaryRegistry bcs;
		std::unique_ptr<topology::TopologicalDOF<HeatEqBundle::NumDOFs>> topoDOF;

		VerificationTopology(Index nx, Index ny, Real x0, Real x1, Real y0, Real y1, Index Px, Index Py) {

			mesh::generator::BlockMesh2D gen{nx, ny, x0, x1, y0, y1, Px, Py};
			mesh = gen.generate();

			auto g = [](Real, const Real*, Real* out) { out[0] = 0.0; };
			using DirichletT = HeatEqBundle::DirichletExpression<decltype(g)>;

			for (Int tag : {0, 1}) {
				auto bc = std::shared_ptr<fem::boundary::BoundaryCondition<DirichletT>>(new fem::boundary::BoundaryCondition<DirichletT>{tag, {fem::boundary::BCCategory::Essential}, DirichletT{g}});
				bcs.registerBC<DirichletT>(bc);
			}

			topoDOF = std::make_unique<topology::TopologicalDOF<HeatEqBundle::NumDOFs>>(mesh, fem::dof::DOFOrdering::Interleaved);

			HeatEqBundle::Basis basis{Px, Py};
			topoDOF->buildConstraints(basis, bcs);

		}

		Index nearestFreeAlgebraicDOF(Real x, Real y) const {

			Index bestNode = 0;
			Real bestDist = std::numeric_limits<Real>::max();

			for (Index nodeID = 0; nodeID < mesh.data.numNodes; ++nodeID) {
				const Real* c = mesh.getNodeCoord(nodeID);
				const Real dist = (c[0] - x) * (c[0] - x) + (c[1] - y) * (c[1] - y);
				if (dist < bestDist) { bestDist = dist; bestNode = nodeID; }
			}

			const Index tdof = topoDOF->getNodeDOF(bestNode, 0);
			return topoDOF->toAlgebraic(tdof);

		}

	}; // struct VerificationTopology

}; // class HeatEquationTransientBackwardEulerTest

TEST_F(HeatEquationTransientBackwardEulerTest, MidpointDecayMatchesAnalyticEigenmodeSolution) {

	auto problem = makeProblem();
	BEStageT stage(*problem);

	solver::config::TimeStepperLoggerConfig loggerCfg;
	loggerCfg.type = solver::config::LoggerConfig::Type::None; // keep test output quiet

	auto stepper = solver::timestepper::makeTimeStepperRunner<BEStageT>(stage, *problem->solverInstance().timestepper, loggerCfg, problem->equationLabel());

	solver::driver::Transient<BEStageT> driver;
	const bool converged = driver.solve(stage, *stepper);

	ASSERT_TRUE(converged);
	EXPECT_NEAR(stepper->currentTime(), tf, 1e-9);

	VerificationTopology verify(nx, ny, x0, x1, y0, y1, Px, Py);
	const Index midDOF = verify.nearestFreeAlgebraicDOF(0.5 * L, 0.5 * (y0 + y1));

	const Real pi = std::acos(Real(-1));
	const Real analytic = std::sin(pi * 0.5 * L / L) * std::exp(-(pi / L) * (pi / L) * tf);
	const Real numeric = problem->solution().data()[midDOF];

	// first-order-in-time (Backward Euler) + Q1-in-space discretization error, at dt=0.01 and a
	// 20-element mesh -- 5% relative tolerance is generous but still catches a broken timestepper
	// wiring (wrong sign, no decay at all, wrong dt, etc.), which is the point of this test.
	EXPECT_NEAR(numeric, analytic, 0.05 * std::abs(analytic));

}

TEST_F(HeatEquationTransientBackwardEulerTest, SolutionDecaysMonotonicallyTowardZero) {

	auto problem = makeProblem();
	BEStageT stage(*problem);

	solver::config::TimeStepperLoggerConfig loggerCfg;
	loggerCfg.type = solver::config::LoggerConfig::Type::None;

	auto stepper = solver::timestepper::makeTimeStepperRunner<BEStageT>(stage, *problem->solverInstance().timestepper, loggerCfg, problem->equationLabel());

	VerificationTopology verify(nx, ny, x0, x1, y0, y1, Px, Py);
	const Index midDOF = verify.nearestFreeAlgebraicDOF(0.5 * L, 0.5 * (y0 + y1));

	Real previous = problem->solution().data()[midDOF];
	EXPECT_GT(previous, 0.0); // sin(pi*0.5) == 1 at the midpoint, t=0

	while (!stepper->finished()) {

		ASSERT_TRUE(stepper->step());

		const Real current = problem->solution().data()[midDOF];
		EXPECT_LT(current, previous); // pure decay, no source, no oscillation expected
		EXPECT_GT(current, 0.0); // never overshoots past zero into negative territory
		previous = current;

	}

}
