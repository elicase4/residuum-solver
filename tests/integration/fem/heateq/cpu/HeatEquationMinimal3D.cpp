#include <gtest/gtest.h>
#include <memory.h>

#include "core/Config.hpp"
#include "core/FEM.hpp"
#include "core/Mesh.hpp"
#include "core/Topology.hpp"
#include "core/Types.hpp"

#include "mesh/ElementFamily.hpp"
#include "mesh/generator/BlockMesh3D.hpp"

#include "equation/heateq/HeatEquation.hpp"

using namespace residuum;

class CPUHeatEquationMinimal3D : public ::testing::Test {
protected:

	// block mesh parameters
	const Real x0 = 0.0;
	const Real x1 = 4.0;
	const Real y0 = 0.0;
	const Real y1 = 4.0;
	const Real z0 = 0.0;
	const Real z1 = 4.0;
	const Index nx = 4;
	const Index ny = 4;
	const Index nz = 4;

	// general problem parameters
	static constexpr Index nsd = 3;
	static constexpr Index Px = 1;
	static constexpr Index Py = 1;
	static constexpr Index Pz = 1;
	static constexpr Index numQuadPoint = 2;

	// specify backend and equation bundle
	using BackendType = linalg::types::backend::CPU;
	using HeatEqBundle = equation::HeatEquation<nsd, 3, mesh::ElementFamily::Hex>;

	// specify basis, quadratures, and element evaluation
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
	fem::boundary::NaturalBoundaryRegistry<HeatEqBundle::EvalQPBdy> NaturalBCRegistry;

	// rhs source function
	static constexpr auto f = [](Real, const Real* x, Real* out){ out[0] = x[0]*x[1]*x[2]; };

	// specify bc functions
	static constexpr auto g = [](Real, const Real* x, Real* out){ out[0] = x[0]; };
	static constexpr auto h = [](Real, const Real* x, Real* out){ out[0] = 0.0; out[1] = 0.0; out[2] = x[0]; };

	// declare assembler
	fem::assembly::Assembler<BackendType> assembler;

	// declare bc applicator
	fem::boundary::BoundaryApplicator<BackendType> bcApplicator;

	// declare models
	HeatEqBundle::DefaultModel defaultModel;
	HeatEqBundle::DefaultModelBdy defaultModelBdy;
	HeatEqBundle::ConductivityModel constantConductivityModel;

	// operator form
	HeatEqBundle::DiffusionForm diffusionForm;
	fem::form::FormRegistry<HeatEqBundle::DiffusionForm> operatorForms{diffusionForm};

	// source form
	HeatEqBundle::SourceFunction<decltype(f)> sourceFunction{f};
	HeatEqBundle::SourceForm<decltype(f)> sourceForm{sourceFunction};
	fem::form::FormRegistry<HeatEqBundle::SourceForm<decltype(f)>> rhsForms{sourceForm};

	// flux form
	HeatEqBundle::FluxBC<decltype(h)> fluxFunction{h};
	HeatEqBundle::FluxForm<decltype(h)> fluxForm{fluxFunction};
	fem::form::FormRegistry<HeatEqBundle::FluxForm<decltype(h)>> naturalBCForms{fluxForm};

	// declare bcs
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc0;
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc1;
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc2;
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc3;
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::FluxBC<decltype(h)>>> bc4;
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::FluxBC<decltype(h)>>> bc5;

	// SetUp method
	void SetUp() override {

		// build 3D block mesh
		mesh3D = gen.generate();

		// create topological DOF manager
		topoDOF3D = std::make_unique<topology::TopologicalDOF<HeatEqBundle::NumDOFs>>(mesh3D, DOFOrdering);

		// set conductivity model parameters
		constantConductivityModel.setConstant(1.0);

		// LEFT/RIGHT/FRONT/BACK essential
		bc0 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{0, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc0);

		bc1 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{1, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc1);

		bc2 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{2, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc2);

		bc3 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{3, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc3);

		// BOTTOM/TOP natural
		bc4 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::FluxBC<decltype(h)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::FluxBC<decltype(h)>>{4, {fem::boundary::BCCategory::Natural}, HeatEqBundle::FluxBC<decltype(h)>{h}});
		NaturalBCRegistry.registerBC<HeatEqBundle::FluxBC<decltype(h)>, decltype(naturalBCForms), HeatEqBundle::DefaultModelBdy>(bc4, naturalBCForms, defaultModelBdy);

		bc5 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::FluxBC<decltype(h)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::FluxBC<decltype(h)>>{5, {fem::boundary::BCCategory::Natural}, HeatEqBundle::FluxBC<decltype(h)>{h}});
		NaturalBCRegistry.registerBC<HeatEqBundle::FluxBC<decltype(h)>, decltype(naturalBCForms), HeatEqBundle::DefaultModelBdy>(bc5, naturalBCForms, defaultModelBdy);

		// build algebraic dofs after all boundaries are registered
		topoDOF3D->buildConstraints(basis, EssentialBCRegistry);

	}

	// global node id for grid index (i,j,k), i,j,k in [0,4]
	Index nodeID(Index i, Index j, Index k) const { return 25*k + 5*j + i; }

};

TEST_F(CPUHeatEquationMinimal3D, DOFHandling){

	EXPECT_EQ(topoDOF3D->numGlobalDOFs(), 125);
	EXPECT_EQ(topoDOF3D->numFreeDOFs(), 45);

	// free: i,j interior (1,2,3), k unrestricted -- BOTTOM/TOP are natural, not essential
	EXPECT_FALSE(topoDOF3D->isConstrained(nodeID(2,2,2)));
	EXPECT_FALSE(topoDOF3D->isConstrained(nodeID(1,1,0)));
	EXPECT_FALSE(topoDOF3D->isConstrained(nodeID(3,3,4)));

	// constrained: touching any lateral (LEFT/RIGHT/FRONT/BACK) face
	EXPECT_TRUE(topoDOF3D->isConstrained(nodeID(0,2,2)));
	EXPECT_TRUE(topoDOF3D->isConstrained(nodeID(4,2,2)));
	EXPECT_TRUE(topoDOF3D->isConstrained(nodeID(2,0,2)));
	EXPECT_TRUE(topoDOF3D->isConstrained(nodeID(2,4,2)));
	EXPECT_TRUE(topoDOF3D->isConstrained(nodeID(0,0,0)));

	// constraint tags match the registered face
	EXPECT_EQ(topoDOF3D->getConstraintTag(nodeID(0,2,2)), 0);
	EXPECT_EQ(topoDOF3D->getConstraintTag(nodeID(4,2,2)), 1);
	EXPECT_EQ(topoDOF3D->getConstraintTag(nodeID(2,0,2)), 2);
	EXPECT_EQ(topoDOF3D->getConstraintTag(nodeID(2,4,2)), 3);

}

TEST_F(CPUHeatEquationMinimal3D, KMatrix){

	Real t = 0.0;

	auto K = assembler.createMatrix<HeatEqBundle::NumDOFs>(mesh3D, *topoDOF3D);
	auto U = assembler.createVector<HeatEqBundle::NumDOFs>(mesh3D, *topoDOF3D);

	EXPECT_EQ(K.nRows(), 45);
	EXPECT_EQ(K.nCols(), 45);
	EXPECT_EQ(U.size(), 45);

	assembler.assembleMatrix<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh3D, *topoDOF3D, t, constantConductivityModel, operatorForms, evalEle, quadVol, U, {nullptr}, K, nullptr);

	const Real tol = 1e-10;

	// diagonal: touching 4 elements at the natural (z) boundary layers, 8 elements interior
	for (Index i = 1; i <= 3; ++i){
		for (Index j = 1; j <= 3; ++j){
			for (Index k = 0; k <= 4; ++k){
				Index a = topoDOF3D->toAlgebraic(nodeID(i,j,k));
				Real expected = (k == 0 || k == 4) ? (4.0/3.0) : (8.0/3.0);
				EXPECT_NEAR(K.data()[K.getDataIndex(a,a)], expected, tol) << "i=" << i << " j=" << j << " k=" << k;
			}
		}
	}

	// off-diagonal spot checks, anchored at the fully-interior node (2,2,2)
	Index center = topoDOF3D->toAlgebraic(nodeID(2,2,2));

	// edge-adjacent (differ in exactly one coordinate): exactly zero for a unit-cube trilinear Hex
	Index edgeNbr = topoDOF3D->toAlgebraic(nodeID(2,2,1));
	EXPECT_NEAR(K.data()[K.getDataIndex(center,edgeNbr)], 0.0, tol);

	// face-diagonal (differ in exactly two coordinates)
	Index faceDiagNbr = topoDOF3D->toAlgebraic(nodeID(1,1,2));
	EXPECT_NEAR(K.data()[K.getDataIndex(center,faceDiagNbr)], -1.0/6.0, tol);

	// body-diagonal (differ in all three coordinates)
	Index bodyDiagNbr = topoDOF3D->toAlgebraic(nodeID(1,1,1));
	EXPECT_NEAR(K.data()[K.getDataIndex(center,bodyDiagNbr)], -1.0/12.0, tol);

}

TEST_F(CPUHeatEquationMinimal3D, OVector){

	Real t = 0.0;

	auto O = assembler.createVector<HeatEqBundle::NumDOFs>(mesh3D, *topoDOF3D);
	auto U = assembler.createVector<HeatEqBundle::NumDOFs>(mesh3D, *topoDOF3D);

	EXPECT_EQ(O.size(), 45);
	EXPECT_EQ(U.size(), 45);

	for (Index i = 0; i < topoDOF3D->numFreeDOFs(); ++i){
		U.data()[i] = 1.0;
	}

	assembler.assembleVector<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh3D, *topoDOF3D, t, constantConductivityModel, operatorForms, evalEle, quadVol, U, nullptr, {nullptr}, O, nullptr);

	const Real tol = 1e-10;

	for (Index k = 0; k <= 4; ++k){
		Index a = topoDOF3D->toAlgebraic(nodeID(2,2,k));
		EXPECT_NEAR(O.data()[a], 0.0, tol) << "k=" << k;
	}

	for (Index i : {1, 3}){
		for (Index j : {1, 3}){
			for (Index k = 0; k <= 4; ++k){
				Index a = topoDOF3D->toAlgebraic(nodeID(i,j,k));
				Real expected = (k == 0 || k == 4) ? (5.0/6.0) : (5.0/3.0);
				EXPECT_NEAR(O.data()[a], expected, tol) << "i=" << i << " j=" << j << " k=" << k;
			}
		}
	}

}

TEST_F(CPUHeatEquationMinimal3D, FVector){

	Real t = 0.0;

	auto F = assembler.createVector<HeatEqBundle::NumDOFs>(mesh3D, *topoDOF3D);
	auto U = assembler.createVector<HeatEqBundle::NumDOFs>(mesh3D, *topoDOF3D);

	EXPECT_EQ(F.size(), 45);
	EXPECT_EQ(U.size(), 45);

	const Real tol = 1e-10;

	// source assembly: f = x*y*z separates as i*j*Sz(k) for unit-cube elements
	const Real Sz[5] = {1.0/6.0, 1.0, 2.0, 3.0, 11.0/6.0};

	assembler.assembleVector<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::DefaultModel, decltype(rhsForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh3D, *topoDOF3D, t, defaultModel, rhsForms, evalEle, quadVol, U, nullptr, {nullptr}, F, nullptr);

	for (Index i = 1; i <= 3; ++i){
		for (Index j = 1; j <= 3; ++j){
			for (Index k = 0; k <= 4; ++k){
				Index a = topoDOF3D->toAlgebraic(nodeID(i,j,k));
				Real expected = static_cast<Real>(i) * static_cast<Real>(j) * Sz[k];
				EXPECT_NEAR(F.data()[a], expected, tol) << "i=" << i << " j=" << j << " k=" << k;
			}
		}
	}

	// apply natural bcs: BOTTOM/TOP flux h=(0,0,x) only reaches the k=0/k=4 layers
	bcApplicator.applyNaturalBCs<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPBdy, HeatEqBundle::QuadratureBoundaryType>(mesh3D, *topoDOF3D, NaturalBCRegistry, t, evalEle, quadBdy, F);

	for (Index i = 1; i <= 3; ++i){
		for (Index j = 1; j <= 3; ++j){
			for (Index k = 0; k <= 4; ++k){
				Index a = topoDOF3D->toAlgebraic(nodeID(i,j,k));
				Real expected = static_cast<Real>(i) * static_cast<Real>(j) * Sz[k];
				if (k == 0) expected -= static_cast<Real>(i);
				if (k == 4) expected += static_cast<Real>(i);
				EXPECT_NEAR(F.data()[a], expected, tol) << "i=" << i << " j=" << j << " k=" << k;
			}
		}
	}

	bcApplicator.applyEssentialBCs<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType>(mesh3D, *topoDOF3D, EssentialBCRegistry, t, constantConductivityModel, operatorForms, evalEle, quadVol, F);

	EXPECT_NEAR(F.data()[topoDOF3D->toAlgebraic(nodeID(1,1,0))], -1.0/3.0, tol);
	EXPECT_NEAR(F.data()[topoDOF3D->toAlgebraic(nodeID(1,1,1))], 2.0, tol);
	EXPECT_NEAR(F.data()[topoDOF3D->toAlgebraic(nodeID(1,1,4))], 10.0/3.0, tol);
	EXPECT_NEAR(F.data()[topoDOF3D->toAlgebraic(nodeID(2,2,2))], 8.0, tol);
	EXPECT_NEAR(F.data()[topoDOF3D->toAlgebraic(nodeID(3,1,1))], 26.0/3.0, tol);

}
