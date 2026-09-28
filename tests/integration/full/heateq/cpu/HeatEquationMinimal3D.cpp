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
#include "mesh/generator/BlockMesh3D.hpp"

#include "equation/heateq/HeatEquation.hpp"

using namespace residuum;

class CPUHeatEquationMinimal3D : public ::testing::Test {
protected:

	// block mesh parameters
	const Real x0 = -1.0;
	const Real x1 = 1.0;
	const Real y0 = -1.0;
	const Real y1 = 1.0;
	const Real z0 = -1.0;
	const Real z1 = 1.0;
	const Index nx = 8;
	const Index ny = 8;
	const Index nz = 8;

	// general mesh parameters
	static constexpr Index nsd = 3;
	static constexpr Index Px = 1;
	static constexpr Index Py = 1;
	static constexpr Index Pz = 1;
	static constexpr Index numQuadPoint = 2;

	// basis, quadrature, and equation type
	using BackendType = linalg::types::backend::CPU;
	using HeatEqBundle = equation::HeatEquation<nsd, 3, mesh::ElementFamily::Hex>;

	// resolved discretization instances
	HeatEqBundle::Basis basis{Px, Py, Pz};
	HeatEqBundle::QuadratureVolumeType quadVol{numQuadPoint, numQuadPoint, numQuadPoint};
	HeatEqBundle::QuadratureBoundaryType quadBdy{numQuadPoint, numQuadPoint};
	HeatEqBundle::EvalEle evalEle{basis};

	// dof parameters
	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Interleaved;

	// initialize mesh and topology
	mesh::generator::BlockMesh3D gen{nx, ny, nz, x0, x1, y0, y1, z0, z1, Px, Py, Pz};
	mesh::Mesh mesh3D;
	std::unique_ptr<topology::TopologicalDOF<HeatEqBundle::NumDOFs>> topoDOF3D;

	// boundary registry
	fem::boundary::EssentialBoundaryRegistry EssentialBCRegistry;

	// rhs source function -- zero, since (1-x)(1-y)(1-z) is harmonic (linear in each variable)
	static constexpr auto f = [](Real, const Real*, Real* out){ out[0] = 0.0; };

	// manufactured solution -- exactly representable by the trilinear Q1 Hex basis
	static constexpr auto g = [](Real, const Real* x, Real* out){ out[0] = (1 - x[0])*(1 - x[1])*(1 - x[2]); };

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

	// declare bcs -- essential on all 6 faces
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc0;
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc1;
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc2;
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc3;
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc4;
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc5;

	// SetUp method
	void SetUp() override {

		// build 3D block mesh
		mesh3D = gen.generate();

		// create topological DOF manager
		topoDOF3D = std::make_unique<topology::TopologicalDOF<HeatEqBundle::NumDOFs>>(mesh3D, DOFOrdering);

		// set conductivity model parameters
		constantConductivityModel.setConstant(1.0);

		// register all 6 faces as essential
		bc0 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{0, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc0);

		bc1 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{1, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc1);

		bc2 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{2, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc2);

		bc3 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{3, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc3);

		bc4 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{4, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc4);

		bc5 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{5, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc5);

		// build algebraic dofs after all boundaries are registered
		topoDOF3D->buildConstraints(basis, EssentialBCRegistry);

	}

};

TEST_F(CPUHeatEquationMinimal3D, DOFHandlingCGSolve){

	// (nx-1)*(ny-1)*(nz-1) interior nodes are free
	EXPECT_EQ(topoDOF3D->numGlobalDOFs(), 729);
	EXPECT_EQ(topoDOF3D->numFreeDOFs(), 343);

}

TEST_F(CPUHeatEquationMinimal3D, MatrixCGSolverTrilinearSolP1){

	Real t = 0.0;

	auto K = assembler.createMatrix<HeatEqBundle::NumDOFs>(mesh3D, *topoDOF3D);
	auto U = assembler.createVector<HeatEqBundle::NumDOFs>(mesh3D, *topoDOF3D);
	auto F = assembler.createVector<HeatEqBundle::NumDOFs>(mesh3D, *topoDOF3D);

	EXPECT_EQ(K.nRows(), 343);
	EXPECT_EQ(K.nCols(), 343);
	EXPECT_EQ(U.size(), 343);
	EXPECT_EQ(F.size(), 343);

	assembler.assembleMatrix<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh3D, *topoDOF3D, t, constantConductivityModel, operatorForms, evalEle, quadVol, U, {nullptr}, K, nullptr);

	assembler.assembleVector<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::DefaultModel, decltype(rhsForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh3D, *topoDOF3D, t, defaultModel, rhsForms, evalEle, quadVol, U, nullptr, {nullptr}, F, nullptr);

	bcApplicator.applyEssentialBCs<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType>(mesh3D, *topoDOF3D, EssentialBCRegistry, t, constantConductivityModel, operatorForms, evalEle, quadVol, F);

	linalg::op::CSROperator<linalg::types::CSRMatrix<Real, BackendType>> op(K);

	linalg::solver::iterative::cg::Workspace<linalg::types::Vector<Real, BackendType>> W(topoDOF3D->numFreeDOFs());
	linalg::solver::SolverReport<linalg::types::Vector<Real, BackendType>> report;

	linalg::solver::preconditioner::Identity<linalg::types::Vector<Real, BackendType>> M;
	utils::logging::linear::ConsoleLogger logger("Heat Equation", "CG", "Identity", {"T"});

	const Real solverTol = 1e-12;
	const Index MaxIter = 10000;
	const linalg::solver::iterative::cg::ToleranceType tolType = linalg::solver::iterative::cg::ToleranceType::Relative;
	linalg::solver::iterative::cg::Config<linalg::types::Vector<Real, BackendType>> cfg{solverTol, tolType, MaxIter};

	linalg::solver::iterative::cg::Solver<decltype(op), decltype(F), decltype(M), decltype(logger)> solver(cfg);

	solver.solve(report, logger, W, M, op, F, U);

	const Real tol = 1e-10;

	EXPECT_TRUE(report.converged);
	EXPECT_NEAR(report.finalResidual, 1e-12, 1e-11);

	Index solIndex = 0;
	for (Index k = 1; k < nz; ++k) {
		for (Index j = 1; j < ny; ++j) {
			for (Index i = 1; i < nx; ++i) {
				Real x_ijk = x0 + i*((x1-x0)/nx);
				Real y_ijk = y0 + j*((y1-y0)/ny);
				Real z_ijk = z0 + k*((z1-z0)/nz);
				EXPECT_NEAR(U.data()[solIndex], (1 - x_ijk)*(1 - y_ijk)*(1 - z_ijk), tol);
				solIndex++;
			}
		}
	}

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "matrixsol3d_output.vtk";
	io::fieldio::FieldIO::writeVTK<HeatEqBundle::NumDOFs>(mesh3D, *topoDOF3D, EssentialBCRegistry, t, U.data(), {"theta"}, path.string());

}

TEST_F(CPUHeatEquationMinimal3D, MatrixFreeCGSolver){

	Real t = 0.0;

	auto U = assembler.createVector<HeatEqBundle::NumDOFs>(mesh3D, *topoDOF3D);
	auto F = assembler.createVector<HeatEqBundle::NumDOFs>(mesh3D, *topoDOF3D);

	EXPECT_EQ(U.size(), 343);
	EXPECT_EQ(F.size(), 343);

	U.zero();

	assembler.assembleVector<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::DefaultModel, decltype(rhsForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh3D, *topoDOF3D, t, defaultModel, rhsForms, evalEle, quadVol, U, nullptr, {nullptr}, F, nullptr);

	bcApplicator.applyEssentialBCs<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType>(mesh3D, *topoDOF3D, EssentialBCRegistry, t, constantConductivityModel, operatorForms, evalEle, quadVol, F);

	linalg::op::FEMOperator<fem::assembly::Assembler<BackendType>, topology::TopologicalDOF<HeatEqBundle::NumDOFs>, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free, linalg::types::Vector<Real, BackendType>> op(assembler, mesh3D, *topoDOF3D, &t, constantConductivityModel, operatorForms, evalEle, quadVol, nullptr, {});

	linalg::solver::iterative::cg::Workspace<linalg::types::Vector<Real, BackendType>> W(topoDOF3D->numFreeDOFs());
	linalg::solver::SolverReport<linalg::types::Vector<Real, BackendType>> report;

	linalg::solver::preconditioner::Identity<linalg::types::Vector<Real, BackendType>> M;
	utils::logging::linear::ConsoleLogger logger("Heat Equation", "CG", "Identity", {"T"});

	const Real solverTol = 1e-12;
	const Index MaxIter = 10000;
	const linalg::solver::iterative::cg::ToleranceType tolType = linalg::solver::iterative::cg::ToleranceType::Relative;
	linalg::solver::iterative::cg::Config<linalg::types::Vector<Real, BackendType>> cfg{solverTol, tolType, MaxIter};

	linalg::solver::iterative::cg::Solver<decltype(op), decltype(F), decltype(M), decltype(logger)> solver(cfg);

	solver.solve(report, logger, W, M, op, F, U);

	const Real tol = 1e-10;

	EXPECT_TRUE(report.converged);
	EXPECT_NEAR(report.finalResidual, 1e-12, 1e-11);

	Index solIndex = 0;
	for (Index k = 1; k < nz; ++k) {
		for (Index j = 1; j < ny; ++j) {
			for (Index i = 1; i < nx; ++i) {
				Real x_ijk = x0 + i*((x1-x0)/nx);
				Real y_ijk = y0 + j*((y1-y0)/ny);
				Real z_ijk = z0 + k*((z1-z0)/nz);
				EXPECT_NEAR(U.data()[solIndex], (1 - x_ijk)*(1 - y_ijk)*(1 - z_ijk), tol);
				solIndex++;
			}
		}
	}

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "matrixfreesol3d_output.vtk";
	io::fieldio::FieldIO::writeVTK<HeatEqBundle::NumDOFs>(mesh3D, *topoDOF3D, EssentialBCRegistry, t, U.data(), {"theta"}, path.string());

}
