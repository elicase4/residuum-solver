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
#include "solver/stage/SteadyStage.hpp"

#include "topology/TopologicalDOF.hpp"

using namespace residuum;

namespace {

	using BackendType = linalg::types::backend::CPU;
	using HeatEqBundle = equation::HeatEquation<2, 2, mesh::ElementFamily::Quad>;
	using HeatProblemT = application::heateq::problem::HeatProblem<BackendType, HeatEqBundle>;
	using SteadyStageT = solver::stage::SteadyStage<HeatProblemT>;

} // namespace

static_assert(solver::stage::Stage<SteadyStageT>);

// SteadyStage carries no equation-specific code of its own (see its header comment) -- these
// tests exist to confirm it correctly orchestrates HeatProblem's toolbox (already checked in
// isolation by UnitTest_Application_Heateq_Problem_HeatProblem) through the generic
// initialize/assemble/solve/finalize lifecycle.
class SteadyStageTest : public ::testing::Test {
protected:

	const Real x0 = -1.0, x1 = 1.0, y0 = -1.0, y1 = 1.0;
	static constexpr Index nx = 4, ny = 4;
	static constexpr Index Px = 1, Py = 1;

	static constexpr Real conductivity = 2.0;
	static constexpr Real a = 1.0; // dT/dx
	static constexpr Real b = 2.0; // dT/dy

	application::heateq::config::HeatConfig makeConfig() const {

		namespace hconfig = application::heateq::config;
		namespace sconfig = solver::config;

		hconfig::HeatConfig cfg;

		cfg.discretization.quadrature.xi = 2;
		cfg.discretization.quadrature.eta = 2;
		cfg.discretization.dofOrdering = fem::dof::DOFOrdering::Interleaved;

		cfg.solver.driver.type = sconfig::DriverConfig::Type::Steady;

		sconfig::LinearSolverConfig linSolver;
		linSolver.type = sconfig::LinearSolverConfig::Type::CG;
		linSolver.operatorType = sconfig::LinearSolverConfig::OperatorType::CSR;
		linSolver.tolerance = 1e-12;
		linSolver.maxIterations = 2000;
		cfg.solver.linear = linSolver;

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

	std::unique_ptr<HeatProblemT> makeProblem() const {

		const auto cfg = makeConfig();

		HeatEqBundle::Basis basis{Px, Py};
		HeatEqBundle::QuadratureVolumeType quadVol{2, 2};
		HeatEqBundle::QuadratureBoundaryType quadBdy{2};

		mesh::generator::BlockMesh2D gen{nx, ny, x0, x1, y0, y1, Px, Py};
		mesh::Mesh mesh = gen.generate();

		return std::make_unique<HeatProblemT>(cfg, std::move(mesh), std::move(basis), std::move(quadVol), std::move(quadBdy));

	}

	// independent of HeatProblem's internal mesh/topoDOF -- rebuilds the same deterministic
	// mesh + constraint layout externally to translate an algebraic DOF index back to a node
	// coordinate for the analytic check below.
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

}; // class SteadyStageTest

TEST_F(SteadyStageTest, LifecycleAssemblesAndSolvesHarmonicDirichletProblemExactly) {

	auto problem = makeProblem();
	SteadyStageT stage(*problem);

	stage.initialize();
	stage.assemble();
	ASSERT_TRUE(stage.solve());
	stage.finalize();

	expectSolutionMatchesAnalyticField(*problem);

}

TEST_F(SteadyStageTest, SolutionAccessorAliasesProblemU) {

	auto problem = makeProblem();
	SteadyStageT stage(*problem);

	stage.assemble();
	ASSERT_TRUE(stage.solve());

	EXPECT_EQ(&stage.solution(), &problem->U());

}
