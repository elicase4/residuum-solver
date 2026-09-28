#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <cmath>
#include <memory.h>

#include "core/Config.hpp"
#include "core/FEM.hpp"
#include "core/IO.hpp"
#include "core/Mesh.hpp"
#include "core/LinAlg.hpp"
#include "core/Topology.hpp"
#include "core/Types.hpp"
#include "core/Utils.hpp"

#include "mesh/ElementFamily.hpp"
#include "mesh/generator/BlockMesh2D.hpp"

#include "equation/heateq/HeatEquation.hpp"

using namespace residuum;

class CPUHeatEquationMinimal : public ::testing::Test {
protected:

	// block mesh parameters
	const Real x0 = -1.0;
	const Real x1 = 1.0;
	const Real y0 = -1.0;
	const Real y1 = 1.0;
	const Index nx = 24;
	const Index ny = 24;

	// general mesh parameters
	static constexpr Index nsd = 2;
	static constexpr Index Px = 1;
	static constexpr Index Py = 1;
	static constexpr Index numQuadPoint = 2;

	// basis, quadrature, and equation type -- Family/NPD are the only
	// compile-time discretization axes now (see the runtime-dispatch
	// refactor); order/quadrature-point-counts are runtime fields on the
	// basis/quadVol/quadBdy instances below.
	using BackendType = linalg::types::backend::CPU;
	using HeatEqBundle = equation::HeatEquation<nsd, 2, mesh::ElementFamily::Quad>;

	// resolved discretization instances -- constructed once, forwarded into
	// every Assembler/BoundaryApplicator/FEMOperator call below (mirrors
	// HeatProblem).
	HeatEqBundle::Basis basis{Px, Py};
	HeatEqBundle::QuadratureVolumeType quadVol{numQuadPoint, numQuadPoint};
	HeatEqBundle::QuadratureBoundaryType quadBdy{numQuadPoint};
	HeatEqBundle::EvalEle evalEle{basis};

	// dof parameters
	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Interleaved;

	// initialize mesh and topology
	mesh::generator::BlockMesh2D gen{nx, ny, x0, x1, y0, y1, Px, Py};
	mesh::Mesh mesh2D;
	std::unique_ptr<topology::TopologicalDOF<HeatEqBundle::NumDOFs>> topoDOF2D;

	// boundary registry
	fem::boundary::EssentialBoundaryRegistry EssentialBCRegistry;

	// rhs source functions
	static constexpr auto f = [](Real, const Real*, Real* out){ out[0] = 0.0; };

	// specify bc functions
	static constexpr auto g = [](Real, const Real* x, Real* out){ out[0] = (1 - x[0])*(1 - x[1]); };

	// declare assembler
	fem::assembly::Assembler<BackendType> assembler;

	// declare bc applicator
	fem::boundary::BoundaryApplicator<BackendType> bcApplicator;

	// declare model
	HeatEqBundle::DefaultModel defaultModel;
	HeatEqBundle::ConductivityModel constantConductivityModel;

	// operator form
	HeatEqBundle::DiffusionForm diffusionForm;
	fem::form::FormRegistry<HeatEqBundle::DiffusionForm> operatorForms{diffusionForm};

	// source form
	HeatEqBundle::SourceFunction<decltype(f)> sourceFunction{f};
	HeatEqBundle::SourceForm<decltype(f)> sourceForm{sourceFunction};
	fem::form::FormRegistry<HeatEqBundle::SourceForm<decltype(f)>> rhsForms{sourceForm};

	// declare bcs
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc0;
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc1;
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc2;
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc3;

	// SetUp method
	void SetUp() override {

		// build 2D block mesh
		mesh2D = gen.generate();

		// create topological DOF manager
		topoDOF2D = std::make_unique<topology::TopologicalDOF<HeatEqBundle::NumDOFs>>(mesh2D, DOFOrdering);

		// set conductivity model parameters
		constantConductivityModel.setConstant(1.0);

		// Set and register boundary 0
		bc0 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{0, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc0);

		// Set and register boundary 1
		bc1 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{1, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc1);

		// Set and register boundary 2
		bc2 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{2, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc2);

		// Set and register boundary 3
		bc3 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{3, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc3);

		// build algrebraic dofs after all boundaries are registered
		topoDOF2D->buildConstraints(basis, EssentialBCRegistry);

	}
};

TEST_F(CPUHeatEquationMinimal, DOFHandlingCGSolve){

	// Test TopologicalDOF
	EXPECT_EQ(topoDOF2D->numGlobalDOFs(), 625);
	EXPECT_EQ(topoDOF2D->numFreeDOFs(), 529);

}

TEST_F(CPUHeatEquationMinimal, MatrixCGSolverBilinearSolP1){

	// arbitrary time
	Real t = 0.0;

	// create system matrix
	auto K = assembler.createMatrix<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF2D);
	auto U = assembler.createVector<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF2D);
	auto F = assembler.createVector<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF2D);

	// test matrix sizes
	EXPECT_EQ(K.nRows(), 529);
	EXPECT_EQ(K.nCols(), 529);

	// test vector sizes
	EXPECT_EQ(U.size(), 529);
	EXPECT_EQ(F.size(), 529);

	// call assembly for system matrix
	assembler.assembleMatrix<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh2D, *topoDOF2D, t, constantConductivityModel, operatorForms, evalEle, quadVol, U, {nullptr}, K, nullptr);

	// call assembly for rhs vector
	assembler.assembleVector<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::DefaultModel, decltype(rhsForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh2D, *topoDOF2D, t, defaultModel, rhsForms, evalEle, quadVol, U, nullptr, {nullptr}, F, nullptr);

	// apply essential bcs
	bcApplicator.applyEssentialBCs<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType>(mesh2D, *topoDOF2D, EssentialBCRegistry, t, constantConductivityModel, operatorForms, evalEle, quadVol, F);

	// define operator
	linalg::op::CSROperator<linalg::types::CSRMatrix<Real, BackendType>> op(K);

	// setup solver workspace & report
	linalg::solver::iterative::cg::Workspace<linalg::types::Vector<Real, BackendType>> W(topoDOF2D->numFreeDOFs());
	linalg::solver::SolverReport<linalg::types::Vector<Real, BackendType>> report;

	// setup preconditioner & logger
	linalg::solver::preconditioner::Identity<linalg::types::Vector<Real, BackendType>> M;
	utils::logging::linear::ConsoleLogger logger("Heat Equation", "CG", "Identity", {"T"});

	// setup solver config
	const Real solverTol = 1e-12;
	const Index MaxIter = 10000;
	const linalg::solver::iterative::cg::ToleranceType tolType = linalg::solver::iterative::cg::ToleranceType::Relative;
	linalg::solver::iterative::cg::Config<linalg::types::Vector<Real, BackendType>> cfg{solverTol, tolType, MaxIter};

	// declare solver
	linalg::solver::iterative::cg::Solver<decltype(op), decltype(F), decltype(M), decltype(logger)> solver(cfg);

	// solver linear system
	solver.solve(report, logger, W, M, op, F, U);

	// test tolerance
	const Real tol = 1e-10;

	// test report
	EXPECT_TRUE(report.converged);
	EXPECT_NEAR(report.finalResidual, 1e-12, 1e-11);

	// test solution
	Index solIndex = 0;
	for (Index j = 1; j < ny; ++j) {
		for (Index i = 1; i < nx; ++i) {
			Real x_ij = x0 + i*((x1-x0)/nx);
			Real y_ij = y0 + j*((y1-y0)/ny);
			EXPECT_NEAR(U.data()[solIndex], (1 - x_ij)*(1 - y_ij), tol);
			solIndex++;
		}
	}

	// write solution
	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "matrixsol_output.vtk";
	io::fieldio::FieldIO::writeVTK<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF2D, EssentialBCRegistry, t, U.data(), {"theta"}, path.string());

}

TEST_F(CPUHeatEquationMinimal, MatrixFreeCGSolver){

	// arbitrary time
	Real t = 0.0;

	// create system matrix
	auto U = assembler.createVector<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF2D);
	auto F = assembler.createVector<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF2D);

	// test vector sizes
	EXPECT_EQ(U.size(), 529);
	EXPECT_EQ(F.size(), 529);

	// fill U
	U.zero();

	// call assembly for rhs vector
	assembler.assembleVector<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::DefaultModel, decltype(rhsForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh2D, *topoDOF2D, t, defaultModel, rhsForms, evalEle, quadVol, U, nullptr, {nullptr}, F, nullptr);

	// apply essential bcs
	bcApplicator.applyEssentialBCs<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType>(mesh2D, *topoDOF2D, EssentialBCRegistry, t, constantConductivityModel, operatorForms, evalEle, quadVol, F);

	// define operator
	linalg::op::FEMOperator<fem::assembly::Assembler<BackendType>, topology::TopologicalDOF<HeatEqBundle::NumDOFs>, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free, linalg::types::Vector<Real, BackendType>> op(assembler, mesh2D, *topoDOF2D, &t, constantConductivityModel, operatorForms, evalEle, quadVol, nullptr, {});

	// setup solver workspace & report
	linalg::solver::iterative::cg::Workspace<linalg::types::Vector<Real, BackendType>> W(topoDOF2D->numFreeDOFs());
	linalg::solver::SolverReport<linalg::types::Vector<Real, BackendType>> report;

	// setup preconditioner & logger
	linalg::solver::preconditioner::Identity<linalg::types::Vector<Real, BackendType>> M;
	utils::logging::linear::ConsoleLogger logger("Heat Equation", "CG", "Identity", {"T"});

	// setup solver config
	const Real solverTol = 1e-12;
	const Index MaxIter = 10000;
	const linalg::solver::iterative::cg::ToleranceType tolType = linalg::solver::iterative::cg::ToleranceType::Relative;
	linalg::solver::iterative::cg::Config<linalg::types::Vector<Real, BackendType>> cfg{solverTol, tolType, MaxIter};

	// declare solver
	linalg::solver::iterative::cg::Solver<decltype(op), decltype(F), decltype(M), decltype(logger)> solver(cfg);

	// solver linear system
	solver.solve(report, logger, W, M, op, F, U);

	// test tolerance
	const Real tol = 1e-10;

	// test report
	EXPECT_TRUE(report.converged);
	EXPECT_NEAR(report.finalResidual, 1e-12, 1e-11);

	// test solution
	Index solIndex = 0;
	for (Index j = 1; j < ny; ++j) {
		for (Index i = 1; i < nx; ++i) {
			Real x_ij = x0 + i*((x1-x0)/nx);
			Real y_ij = y0 + j*((y1-y0)/ny);
			EXPECT_NEAR(U.data()[solIndex], (1 - x_ij)*(1 - y_ij), tol);
			solIndex++;
		}
	}

	// write solution
	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "matrixfreesol_output.vtk";
	io::fieldio::FieldIO::writeVTK<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF2D, EssentialBCRegistry, t, U.data(), {"theta"}, path.string());

}
