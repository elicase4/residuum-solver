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

#include "solver/stage/SteadyStage.hpp"

#include "topology/TopologicalDOF.hpp"

using namespace residuum;

//   u(x)     = u(theta0) + [u(theta1) - u(theta0)] * (x - x0)/(x1 - x0)
//   theta(x) = [-k0 + sqrt(k0^2 + 2*k1*u(x))] / k1

class HeatEquationSteadyNonlinearConductivityTest : public ::testing::Test {
protected:

	const Real x0 = 0.0, x1 = 1.0, y0 = 0.0, y1 = 1.0;
	static constexpr Index nx = 40, ny = 4;
	static constexpr Index Px = 1, Py = 1;

	static constexpr Real k0 = 1.0;
	static constexpr Real k1 = 0.5;
	static constexpr Real theta0 = 0.0;
	static constexpr Real theta1 = 1.0;

	using BackendType = linalg::types::backend::CPU;
	using HeatEqBundle = equation::HeatEquation<2, 2, mesh::ElementFamily::Quad>;
	using HeatProblemT = application::heateq::problem::HeatProblem<BackendType, HeatEqBundle>;
	using SteadyStageT = solver::stage::SteadyStage<HeatProblemT>;

	static Real kirchhoff(Real theta) {
		return k0 * theta + 0.5 * k1 * theta * theta;
	}

	Real analyticTheta(Real x) const {
		const Real s = (x - x0) / (x1 - x0);
		const Real u = kirchhoff(theta0) + (kirchhoff(theta1) - kirchhoff(theta0)) * s;
		return (-k0 + std::sqrt(k0 * k0 + 2.0 * k1 * u)) / k1;
	}

	application::heateq::config::HeatConfig makeConfig(solver::config::LinearSolverConfig::OperatorType opType = solver::config::LinearSolverConfig::OperatorType::CSR) const {

		namespace hconfig = application::heateq::config;
		namespace sconfig = solver::config;

		hconfig::HeatConfig cfg;

		cfg.discretization.quadrature.xi = 2;
		cfg.discretization.quadrature.eta = 2;
		cfg.discretization.dofOrdering = fem::dof::DOFOrdering::Interleaved;

		cfg.solver.driver.type = sconfig::DriverConfig::Type::Steady;

		sconfig::NonlinearSolverConfig nlSolver;
		nlSolver.type = sconfig::NonlinearSolverConfig::Type::Newton;
		nlSolver.absoluteTolerance = 1e-10;
		nlSolver.relativeTolerance = 1e-10;
		nlSolver.maxIterations = 50;
		nlSolver.linearSolver.type = sconfig::LinearSolverConfig::Type::CG;
		nlSolver.linearSolver.operatorType = opType;
		nlSolver.linearSolver.tolerance = 1e-12;
		nlSolver.linearSolver.maxIterations = 5000;
		cfg.solver.nonlinear = nlSolver;

		cfg.conductivity.type = hconfig::ConductivityConfig::Type::TemperatureDependentIsotropic;
		cfg.conductivity.valueExpression = "1.0 + 0.5*T";
		cfg.conductivity.gradientExpression = "0.5";
		cfg.conductivity.unit = "W/(m*K)";

		cfg.source.type = hconfig::SourceConfig::Type::VolumetricHeatSource;
		cfg.source.read.mode = sconfig::NodalFieldReadConfig::Mode::Expression;
		cfg.source.read.expression = "0.0";
		cfg.source.read.unit = "W/m^3";

		cfg.initialCondition.read.mode = sconfig::NodalFieldReadConfig::Mode::Expression;
		cfg.initialCondition.read.expression = "0.0";
		cfg.initialCondition.read.unit = "K";

		for (Int tag : {0, 1}) {
			hconfig::BoundaryConditionConfig bc;
			bc.boundaryID = tag;
			bc.type = hconfig::BoundaryConditionConfig::Type::Value;
			bc.mode = sconfig::NodalFieldReadConfig::Mode::Expression;
			bc.expression = (tag == 0) ? "0.0" : "1.0";
			bc.unit = "K";
			bc.forms = {hconfig::BoundaryConditionConfig::Form::ValueBC};
			bc.model = hconfig::ConductivityConfig::Type::TemperatureDependentIsotropic;
			cfg.boundaryConditions.push_back(bc);
		}

		return cfg;

	}

	std::unique_ptr<HeatProblemT> makeProblem(solver::config::LinearSolverConfig::OperatorType opType = solver::config::LinearSolverConfig::OperatorType::CSR) const {

		const auto cfg = makeConfig(opType);

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

	void verifyMatchesAnalytic(const HeatProblemT& problem) const {

		VerificationTopology verify(nx, ny, x0, x1, y0, y1, Px, Py);

		const Real tol = 1e-6;
		for (Index i = 1; i < nx; ++i) {
			const Real x_i = x0 + i * ((x1 - x0) / nx);
			const Real expected = analyticTheta(x_i);
			for (Real yFrac : {0.0, 0.5, 1.0}) {
				const Real y_j = y0 + yFrac * (y1 - y0);
				const Index dof = verify.nearestFreeAlgebraicDOF(x_i, y_j);
				EXPECT_NEAR(problem.solution().data()[dof], expected, tol) << "x=" << x_i << " y=" << y_j;
			}
		}

	}

}; // class HeatEquationSteadyNonlinearConductivityTest

TEST_F(HeatEquationSteadyNonlinearConductivityTest, ConvergesToKirchhoffClosedFormSolution) {

	auto problem = makeProblem();
	SteadyStageT stage(*problem);

	stage.assemble();
	ASSERT_TRUE(stage.solve());

	verifyMatchesAnalytic(*problem);

}

// real lifted field for K(theta)/dK/dtheta evaluation) that the CSR path never touches.
TEST_F(HeatEquationSteadyNonlinearConductivityTest, ConvergesWithMatrixFreeFEMOperator) {

	auto problem = makeProblem(solver::config::LinearSolverConfig::OperatorType::FEM);
	SteadyStageT stage(*problem);

	stage.assemble();
	ASSERT_TRUE(stage.solve());

	verifyMatchesAnalytic(*problem);

}
