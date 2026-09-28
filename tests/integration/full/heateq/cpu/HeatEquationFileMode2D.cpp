#include <filesystem>
#include <gtest/gtest.h>

#include "application/heateq/config/HeatConfig.hpp"
#include "application/heateq/problem/HeatProblem.hpp"

#include "core/Types.hpp"
#include "io/fieldio/FieldIO.hpp"

#include "equation/heateq/HeatEquation.hpp"
#include "mesh/ElementFamily.hpp"
#include "mesh/generator/BlockMesh2D.hpp"
#include "solver/stage/SteadyStage.hpp"

using namespace residuum;

class CPUHeatEquationFileMode : public ::testing::Test {
protected:

	const Real x0 = -1.0;
	const Real x1 = 1.0;
	const Real y0 = -1.0;
	const Real y1 = 1.0;
	const Index nx = 24;
	const Index ny = 24;

	static constexpr Index nsd = 2;
	static constexpr Index Px = 1;
	static constexpr Index Py = 1;
	static constexpr Index numQuadPoint = 2;

	using BackendType = linalg::types::backend::CPU;
	using HeatEqBundle = equation::HeatEquation<nsd, 2, mesh::ElementFamily::Quad>;

	mesh::generator::BlockMesh2D gen{nx, ny, x0, x1, y0, y1, Px, Py};

};

TEST_F(CPUHeatEquationFileMode, AllDirichletFileModeMatchesKnownSolution) {

	mesh::Mesh mesh = gen.generate();

	std::vector<Real> tField(mesh.data.numNodes * HeatEqBundle::NumDOFs);
	std::vector<Real> zeroField(mesh.data.numNodes * HeatEqBundle::NumDOFs, Real(0));
	for (Index nodeID = 0; nodeID < mesh.data.numNodes; ++nodeID) {
		const Real* c = mesh.getNodeCoord(nodeID);
		tField[nodeID] = (1 - c[0]) * (1 - c[1]);
	}

	const auto tPath = (std::filesystem::path(TEST_OUTPUT_PATH) / "filemode_dirichlet.pndf").string();
	const auto zeroPath = (std::filesystem::path(TEST_OUTPUT_PATH) / "filemode_zero.pndf").string();
	io::fieldio::FieldIO::writeBinaryRaw<HeatEqBundle::NumDOFs>(mesh, 0.0, tField, tPath);
	io::fieldio::FieldIO::writeBinaryRaw<HeatEqBundle::NumDOFs>(mesh, 0.0, zeroField, zeroPath);

	application::heateq::config::HeatConfig cfg;

	cfg.conductivity.type = application::heateq::config::ConductivityConfig::Type::Constant;
	cfg.conductivity.value = 1.0;

	cfg.source.type = application::heateq::config::SourceConfig::Type::VolumetricHeatSource;
	cfg.source.read.mode = solver::config::NodalFieldReadConfig::Mode::File;
	cfg.source.read.file = zeroPath;

	cfg.initialCondition.read.mode = solver::config::NodalFieldReadConfig::Mode::File;
	cfg.initialCondition.read.file = zeroPath;

	for (Int boundaryID = 0; boundaryID < 4; ++boundaryID) {

		application::heateq::config::BoundaryConditionConfig bc;
		bc.boundaryID = boundaryID;
		bc.type = application::heateq::config::BoundaryConditionConfig::Type::Value;
		bc.mode = solver::config::NodalFieldReadConfig::Mode::File;
		bc.file = tPath;
		bc.forms = {application::heateq::config::BoundaryConditionConfig::Form::ValueBC};
		bc.model = application::heateq::config::ConductivityConfig::Type::Constant;
		cfg.boundaryConditions.push_back(bc);

	}

	cfg.solver.driver.type = solver::config::DriverConfig::Type::Steady;
	cfg.solver.linear = solver::config::LinearSolverConfig{};
	cfg.solver.linear->operatorType = solver::config::LinearSolverConfig::OperatorType::CSR;
	cfg.solver.linear->tolerance = 1e-12;
	cfg.solver.linear->maxIterations = 10000;

	HeatEqBundle::Basis basis{Px, Py};
	HeatEqBundle::QuadratureVolumeType quadVol{numQuadPoint, numQuadPoint};
	HeatEqBundle::QuadratureBoundaryType quadBdy{numQuadPoint};

	application::heateq::problem::HeatProblem<BackendType, HeatEqBundle> problem(cfg, std::move(mesh), basis, quadVol, quadBdy);

	EXPECT_EQ(problem.numFreeDOFs(), 529);

	solver::stage::SteadyStage<application::heateq::problem::HeatProblem<BackendType, HeatEqBundle>> stage(problem);
	stage.assemble();
	EXPECT_TRUE(stage.solve());

	const Real tol = 1e-9;
	Index solIndex = 0;
	for (Index j = 1; j < ny; ++j) {
		for (Index i = 1; i < nx; ++i) {
			Real x_ij = x0 + i*((x1-x0)/nx);
			Real y_ij = y0 + j*((y1-y0)/ny);
			EXPECT_NEAR(problem.solution().data()[solIndex], (1 - x_ij)*(1 - y_ij), tol);
			solIndex++;
		}
	}

}
