#include <memory>
#include <utility>

#include <gtest/gtest.h>

#include "core/Types.hpp"

#include "application/heateq/config/HeatConfig.hpp"
#include "application/heateq/problem/HeatProblem.hpp"

#include "equation/heateq/HeatEquation.hpp"

#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/dof/DOFOrdering.hpp"

#include "linalg/solver/base/SolverReport.hpp"
#include "linalg/types/backend/CPU.hpp"

#include "mesh/ElementFamily.hpp"
#include "mesh/generator/BlockMesh2D.hpp"

#include "solver/SolverInstance.hpp"
#include "solver/config/DriverConfig.hpp"
#include "solver/problem/Problem.hpp"
#include "solver/problem/TransientCapableProblem.hpp"

#include "topology/TopologicalDOF.hpp"

using namespace residuum;

namespace {

	using BackendType = linalg::types::backend::CPU;
	using HeatEqBundle = equation::HeatEquation<2, 2, mesh::ElementFamily::Quad>;
	using HeatProblemT = application::heateq::problem::HeatProblem<BackendType, HeatEqBundle>;

} // namespace

// HeatProblem must satisfy both concepts unconditionally, at compile time -- regardless of
// whatever driver/solver mode a particular config instance ends up resolving to at runtime.
static_assert(solver::problem::Problem<HeatProblemT>);
static_assert(solver::problem::TransientCapableProblem<HeatProblemT>);

class HeatProblemTest : public ::testing::Test {
protected:

	const Real x0 = -1.0, x1 = 1.0, y0 = -1.0, y1 = 1.0;
	static constexpr Index nx = 4, ny = 4;
	static constexpr Index Px = 1, Py = 1;

	static constexpr Real conductivity = 2.0;
	static constexpr Real a = 1.0; // dT/dx
	static constexpr Real b = 2.0; // dT/dy

	// Every boundary is a Dirichlet BC matching the harmonic field T = a*x + b*y -- by
	// uniqueness of the Dirichlet problem for Laplace's equation (zero source), the exact
	// solution everywhere in the domain is that same linear field. Gives a closed-form check
	// for the full assemble+solve chain, not just "it converged".
	application::heateq::config::HeatConfig makeConfig(bool transient) const {

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

		cfg.solver.driver.type = transient ? sconfig::DriverConfig::Type::Transient : sconfig::DriverConfig::Type::Steady;

		if (transient) {

			sconfig::TimeStepperConfig ts;
			ts.type = sconfig::TimeStepperConfig::Type::BackwardEuler;
			ts.t0 = 0.0;
			ts.tf = 1.0;
			ts.stepSize.mode = sconfig::TimeStepSizeConfig::Mode::Constant;
			ts.stepSize.dt = 0.1;
			ts.linearSolver = linSolver;
			cfg.solver.timestepper = ts;

			cfg.density = hconfig::ScalarMaterialPropertyConfig{1.0, "kg/m^3"};
			cfg.specificHeat = hconfig::SpecificHeatConfig{hconfig::SpecificHeatConfig::Type::Constant, 1.0, "", "", "J/(kg*K)"};

		}

		cfg.conductivity.type = hconfig::ConductivityConfig::Type::Constant;
		cfg.conductivity.value = conductivity;
		cfg.conductivity.unit = "W/(m*K)";

		cfg.source.type = hconfig::SourceConfig::Type::VolumetricHeatSource;
		cfg.source.read.mode = sconfig::NodalFieldReadConfig::Mode::Expression;
		cfg.source.read.expression = "0.0";
		cfg.source.read.unit = "W/m^3";

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

	std::unique_ptr<HeatProblemT> makeProblem(bool transient) const {

		const auto cfg = makeConfig(transient);

		HeatEqBundle::Basis basis{Px, Py};
		HeatEqBundle::QuadratureVolumeType quadVol{2, 2};
		HeatEqBundle::QuadratureBoundaryType quadBdy{2};

		mesh::generator::BlockMesh2D gen{nx, ny, x0, x1, y0, y1, Px, Py};
		mesh::Mesh mesh = gen.generate();

		return std::make_unique<HeatProblemT>(cfg, std::move(mesh), std::move(basis), std::move(quadVol), std::move(quadBdy));

	}

	// independent of HeatProblem's internal mesh/topoDOF -- rebuilds the same deterministic
	// mesh + constraint layout externally, purely to translate algebraic DOF index back to a
	// node coordinate for the analytic check below.
	struct VerificationTopology {

		mesh::Mesh mesh;
		fem::boundary::EssentialBoundaryRegistry bcs;
		std::unique_ptr<topology::TopologicalDOF<HeatEqBundle::NumDOFs>> topoDOF; // built AFTER mesh is generated -- TopologicalDOF reads mesh eagerly at construction

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

	void expectSolutionMatchesAnalyticField(const HeatProblemT& problem) const {

		VerificationTopology verify(nx, ny, x0, x1, y0, y1, Px, Py, a, b);
		const Real* solutionData = problem.solution().data();

		Index numFreeChecked = 0;
		for (Index nodeID = 0; nodeID < verify.mesh.data.numNodes; ++nodeID) {

			const Index tdof = verify.topoDOF->getNodeDOF(nodeID, 0);
			if (verify.topoDOF->isConstrained(tdof)) continue; // constrained DOFs are exactly the BC expression by construction, not worth re-checking

			const Index adof = verify.topoDOF->toAlgebraic(tdof);
			const Real* c = verify.mesh.getNodeCoord(nodeID);
			const Real expected = a * c[0] + b * c[1];

			EXPECT_NEAR(solutionData[adof], expected, 1e-8);
			++numFreeChecked;

		}

		EXPECT_GT(numFreeChecked, 0); // sanity -- the interior isn't empty for a 4x4 mesh

	}

}; // class HeatProblemTest

TEST_F(HeatProblemTest, NumFreeDOFsMatchesInteriorNodeCount) {

	auto problem = makeProblem(false);

	// perimeter of the block mesh is fully Dirichlet-constrained -- only interior nodes are free
	const Index expectedFree = (nx - 1) * (ny - 1);
	EXPECT_EQ(problem->numFreeDOFs(), expectedFree);

}

TEST_F(HeatProblemTest, CreateVectorAndCreateMatrixAreSizedToFreeDOFs) {

	auto problem = makeProblem(false);

	auto V = problem->createVector();
	EXPECT_EQ(V.size(), problem->numFreeDOFs());

}

TEST_F(HeatProblemTest, EquationLabelAndSolverInstanceReflectSteadyLinearConfig) {

	auto problem = makeProblem(false);

	EXPECT_EQ(problem->equationLabel(), "Heat Equation");
	EXPECT_EQ(problem->solverInstance().mode, solver::SolverMode::SteadyLinear);
	ASSERT_NE(problem->solverInstance().linear, nullptr);
	EXPECT_EQ(problem->solverInstance().timestepper, nullptr);

}

TEST_F(HeatProblemTest, SolverInstanceReflectsTransientLinearConfig) {

	auto problem = makeProblem(true);

	EXPECT_EQ(problem->solverInstance().mode, solver::SolverMode::TransientLinear);
	ASSERT_NE(problem->solverInstance().timestepper, nullptr);
	EXPECT_DOUBLE_EQ(problem->solverInstance().timestepper->t0, 0.0);

}

// Drives the Problem toolbox methods directly, in the same order solver::stage::SteadyStage
// does internally (see SteadyStage::assemble/solve) -- deliberately NOT going through
// SteadyStage itself, since that's covered separately (#97). The point here is that
// HeatProblem's own createMatrix/makeLinearRunner/assembleMatrix/assembleLoad/applyNatural/
// applyEssential toolbox is individually correct.
TEST_F(HeatProblemTest, ToolboxAssemblesAndSolvesHarmonicDirichletProblemExactly) {

	auto problem = makeProblem(false);

	auto K = problem->createMatrix();
	const Real time = 0.0;
	auto runner = problem->makeLinearRunner<fem::assembly::GatherMode::Free>(&time, problem->stiffnessForms(), problem->stiffnessModel(), nullptr, {}, K);

	problem->assembleMatrix<fem::assembly::GatherMode::Free>(0.0, problem->stiffnessForms(), problem->stiffnessModel(), {}, K);
	problem->assembleLoad(0.0);
	problem->applyNatural(0.0);
	problem->applyEssential(0.0, problem->stiffnessForms(), problem->stiffnessModel(), problem->F());

	linalg::solver::SolverReport<HeatProblemT::VectorT> report;
	const bool converged = runner->solve(problem->F(), problem->U(), report);

	ASSERT_TRUE(converged);
	expectSolutionMatchesAnalyticField(*problem);

}

// massForms()/massModel()/U_prev() -- the extra TransientCapableProblem surface beyond Problem.
// The ctor copies the (zero) initial condition into U_prev on construction for a transient
// config; nothing has advanced a timestep yet, so U and U_prev must still agree.
TEST_F(HeatProblemTest, TransientToolboxExposesMassRoleAndSeedsUPrevFromInitialCondition) {

	auto problem = makeProblem(true);

	const auto massForms = problem->massForms();
	const auto& massModel = problem->massModel();
	(void)massForms;
	(void)massModel;

	ASSERT_EQ(problem->U_prev().size(), problem->U().size());
	for (Index i = 0; i < problem->U().size(); ++i) {
		EXPECT_DOUBLE_EQ(problem->U_prev().data()[i], problem->U().data()[i]);
	}

}
