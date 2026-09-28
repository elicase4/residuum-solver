#include <gtest/gtest.h>
#include <memory.h>

#include "core/Config.hpp"
#include "core/FEM.hpp"
#include "core/Mesh.hpp"
#include "core/Topology.hpp"
#include "core/Types.hpp"

#include "mesh/ElementFamily.hpp"
#include "mesh/generator/BlockMesh2D.hpp"

#include "equation/heateq/HeatEquation.hpp"

using namespace residuum;

class CPUHeatEquationMinimal : public ::testing::Test {
protected:

	// block mesh parameters
	const Real x0 = 0.0;
	const Real x1 = 4.0;
	const Real y0 = 0.0;
	const Real y1 = 4.0;
	const Index nx = 4;
	const Index ny = 4;

	// general problem parameters
	static constexpr Index nsd = 2;
	static constexpr Index Px = 1;
	static constexpr Index Py = 1;
	static constexpr Index numQuadPoint = 2;

	// specify backend and equation bundle
	using BackendType = linalg::types::backend::CPU;
	using HeatEqBundle = equation::HeatEquation<nsd, 2, mesh::ElementFamily::Quad>;

	// specify basis, quadratures, and element evalaution
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
	fem::boundary::NaturalBoundaryRegistry<HeatEqBundle::EvalQPBdy> NaturalBCRegistry;

	// rhs source functions
	static constexpr auto f = [](Real, const Real* x, Real* out){ out[0] = x[0]*x[1]; };

	// specify bc functions
	static constexpr auto g = [](Real, const Real* x, Real* out){ out[0] = x[0]; };
	static constexpr auto h = [](Real, const Real* x, Real* out){ out[0] = 0.0; out[1] = x[0]; };

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
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::FluxBC<decltype(h)>>> bc2;
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
		bc2 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::FluxBC<decltype(h)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::FluxBC<decltype(h)>>{2, {fem::boundary::BCCategory::Natural}, HeatEqBundle::FluxBC<decltype(h)>{h}});
		NaturalBCRegistry.registerBC<HeatEqBundle::FluxBC<decltype(h)>, decltype(naturalBCForms), HeatEqBundle::DefaultModelBdy>(bc2, naturalBCForms, defaultModelBdy);

		// Set and register boundary 3
		bc3 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{3, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc3);

		// build algrebraic dofs after all boundaries are registered
		topoDOF2D->buildConstraints(basis, EssentialBCRegistry);

	}
};

TEST_F(CPUHeatEquationMinimal, DOFHandling){

	// Test topologicalDOF
	EXPECT_EQ(topoDOF2D->numGlobalDOFs(), 25);
	EXPECT_EQ(topoDOF2D->numFreeDOFs(), 12);

	// Test topological DOF to algebraic DOF mapping
	for (Index i = 0; i <= 20; i+=5){
		EXPECT_EQ(topoDOF2D->toAlgebraic(i), -1);
	}
	for (Index i = 4; i <= 24; i+=5){
		EXPECT_EQ(topoDOF2D->toAlgebraic(i), -1);
	}
	for (Index i = 21; i <= 22; i++){
		EXPECT_EQ(topoDOF2D->toAlgebraic(i), -1);
	}
	EXPECT_EQ(topoDOF2D->toAlgebraic(1),  0);
	EXPECT_EQ(topoDOF2D->toAlgebraic(2),  1);
	EXPECT_EQ(topoDOF2D->toAlgebraic(3),  2);
	EXPECT_EQ(topoDOF2D->toAlgebraic(6),  3);
	EXPECT_EQ(topoDOF2D->toAlgebraic(7),  4);
	EXPECT_EQ(topoDOF2D->toAlgebraic(8),  5);
	EXPECT_EQ(topoDOF2D->toAlgebraic(11), 6);
	EXPECT_EQ(topoDOF2D->toAlgebraic(12), 7);
	EXPECT_EQ(topoDOF2D->toAlgebraic(13), 8);
	EXPECT_EQ(topoDOF2D->toAlgebraic(16), 9);
	EXPECT_EQ(topoDOF2D->toAlgebraic(17), 10);
	EXPECT_EQ(topoDOF2D->toAlgebraic(18), 11);

	// Test algebraic DOF to topological DOF mapping
	EXPECT_EQ(topoDOF2D->toTopological(0),  1);
	EXPECT_EQ(topoDOF2D->toTopological(1),  2);
	EXPECT_EQ(topoDOF2D->toTopological(2),  3);
	EXPECT_EQ(topoDOF2D->toTopological(3),  6);
	EXPECT_EQ(topoDOF2D->toTopological(4),  7);
	EXPECT_EQ(topoDOF2D->toTopological(5),  8);
	EXPECT_EQ(topoDOF2D->toTopological(6),  11);
	EXPECT_EQ(topoDOF2D->toTopological(7),  12);
	EXPECT_EQ(topoDOF2D->toTopological(8),  13);
	EXPECT_EQ(topoDOF2D->toTopological(9),  16);
	EXPECT_EQ(topoDOF2D->toTopological(10), 17);
	EXPECT_EQ(topoDOF2D->toTopological(11), 18);

	// Test constraint flag
	EXPECT_TRUE(topoDOF2D->isConstrained(15));
	EXPECT_TRUE(topoDOF2D->isConstrained(22));
	EXPECT_TRUE(topoDOF2D->isConstrained(14));
	EXPECT_FALSE(topoDOF2D->isConstrained(17));
	EXPECT_FALSE(topoDOF2D->isConstrained(6));

	// Test get constraint tag
	EXPECT_EQ(topoDOF2D->getConstraintTag(10), 0);
	EXPECT_EQ(topoDOF2D->getConstraintTag(15), 0);
	EXPECT_EQ(topoDOF2D->getConstraintTag(21), 3);
	EXPECT_EQ(topoDOF2D->getConstraintTag(9), 1);

}

TEST_F(CPUHeatEquationMinimal, KMatrix){

	// arbitrary time
	Real t = 0.0;

	// create system matrix
	auto K = assembler.createMatrix<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF2D);
	auto U = assembler.createVector<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF2D);

	// test matrix sizes
	EXPECT_EQ(K.nRows(), 12);
	EXPECT_EQ(K.nCols(), 12);

	// test vector sizes
	EXPECT_EQ(U.size(), 12);

	// call assembly for system matrix
	assembler.assembleMatrix<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh2D, *topoDOF2D, t, constantConductivityModel, operatorForms, evalEle, quadVol, U, {nullptr}, K, nullptr);

	// test tolerance
	const Real tol = 1e-10;

	// tests for system matrix

	// diagonal: bottom-free nodes (touching 2 elements)
	for (Index i = 0; i < 3; ++i){
		EXPECT_NEAR(K.data()[K.getDataIndex(i,i)], 4.0/3.0, tol);
	}

	// diagonal: interior-free nodes (touching 2 elements)
	for (Index i = 3; i < 12; ++i){
		EXPECT_NEAR(K.data()[K.getDataIndex(i,i)], 8.0/3.0, tol);
	}

	// off-diagonal: same-row adjacent (share one element edge - x-adjacent with no upper shared element)
	for (Index i = 0; i < 2; ++i){
		EXPECT_NEAR(K.data()[K.getDataIndex(i,i+1)], -1.0/6.0, tol);
		EXPECT_NEAR(K.data()[K.getDataIndex(i+1,i)], -1.0/6.0, tol);
	}

	// off-diagonal: interior-free nodes (nodes sharing 2 elements - x adjacent and y-adjacent interior pairs)
	for (Index i = 3; i < 4; ++i){
		EXPECT_NEAR(K.data()[K.getDataIndex(i,i+1)], -1.0/3.0, tol);
		EXPECT_NEAR(K.data()[K.getDataIndex(i+1,i)], -1.0/3.0, tol);
	}
	for (Index i = 6; i < 8; ++i){
		EXPECT_NEAR(K.data()[K.getDataIndex(i,i+1)], -1.0/3.0, tol);
		EXPECT_NEAR(K.data()[K.getDataIndex(i+1,i)], -1.0/3.0, tol);
	}
	for (Index i = 9; i < 11; ++i){
		EXPECT_NEAR(K.data()[K.getDataIndex(i,i+1)], -1.0/3.0, tol);
		EXPECT_NEAR(K.data()[K.getDataIndex(i+1,i)], -1.0/3.0, tol);
	}

	// off-diagonal: same column-adjacent (nodes sharing 2 elements - y-adjacent pairs)
	for (Index i = 0; i < 9; ++i){
		EXPECT_NEAR(K.data()[K.getDataIndex(i,i+3)], -1.0/3.0, tol);
		EXPECT_NEAR(K.data()[K.getDataIndex(i+3,i)], -1.0/3.0, tol);
	}

}

TEST_F(CPUHeatEquationMinimal, OVector){

	// arbitrary time
	Real t = 0.0;

	auto O = assembler.createVector<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF2D);
	auto U = assembler.createVector<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF2D);

	// test vector sizes
	EXPECT_EQ(O.size(), 12);
	EXPECT_EQ(U.size(), 12);

	// fill U
	for (Index i = 0; i < topoDOF2D->numFreeDOFs(); ++i){
		U.data()[i] = 1.0;
	}

	// call assembly for system matrix
	assembler.assembleVector<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh2D, *topoDOF2D, t, constantConductivityModel, operatorForms, evalEle, quadVol, U, nullptr, {nullptr}, O, nullptr);

	// test tolerance
	const Real tol = 1e-10;

	// tests for system operator

	// bottom-row corner nodes
	EXPECT_NEAR(O.data()[0], 0.5,     tol); // row_sum = 4/3 - 1/6 - 1/3 - 1/3 = 1/2
	EXPECT_NEAR(O.data()[2], 0.5,     tol); // row_sum = 4/3 - 1/6 - 1/3 - 1/3 = 1/2

	// bottom-row mid node
	EXPECT_NEAR(O.data()[1], 0.0,     tol); // row_sum = 0

	// left/right column interior nodes
	EXPECT_NEAR(O.data()[3], 1.0,     tol); // row_sum = 1
	EXPECT_NEAR(O.data()[5], 1.0,     tol); // row_sum = 1
	EXPECT_NEAR(O.data()[6], 1.0,     tol); // row_sum = 1
	EXPECT_NEAR(O.data()[8], 1.0,     tol); // row_sum = 1

	// full interior nodes
	EXPECT_NEAR(O.data()[4], 0.0,     tol); // row_sum = 0
	EXPECT_NEAR(O.data()[7], 0.0,     tol); // row_sum = 0

	// bottom-row corner nodes
	EXPECT_NEAR(O.data()[9], 5.0/3.0,  tol); // row_sum = 8/3 - 1/3 - 1/3 - 1/3 = 5/3
	EXPECT_NEAR(O.data()[11], 5.0/3.0,  tol); // row_sum = 8/3 - 1/3 - 1/3 - 1/3 = 5/3

	// top-row mid node
	EXPECT_NEAR(O.data()[10], 1.0,     tol); // row_sum = 8/3 - 1/3 - 1/3 = 1

}

TEST_F(CPUHeatEquationMinimal, FVector){

	// arbitrary time
	Real t = 0.0;

	// create system matrix
	auto F = assembler.createVector<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF2D);
	auto U = assembler.createVector<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF2D);

	// test vector sizes
	EXPECT_EQ(F.size(), 12);
	EXPECT_EQ(U.size(), 12);

	// test tolerance
	const Real tol = 1e-10;

	// call assembly for force vector
	assembler.assembleVector<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::DefaultModel, decltype(rhsForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh2D, *topoDOF2D, t, defaultModel, rhsForms, evalEle, quadVol, U, nullptr, {nullptr}, F, nullptr);

	// test before bc application
	EXPECT_NEAR(F.data()[0], 1.0/6.0, tol);
	EXPECT_NEAR(F.data()[1], 1.0/3.0, tol);
	EXPECT_NEAR(F.data()[2], 1.0/2.0, tol);
	EXPECT_NEAR(F.data()[3], 1.0, tol);
	EXPECT_NEAR(F.data()[4], 2.0, tol);
	EXPECT_NEAR(F.data()[5], 3.0, tol);
	EXPECT_NEAR(F.data()[6], 2.0, tol);
	EXPECT_NEAR(F.data()[7], 4.0, tol);
	EXPECT_NEAR(F.data()[8], 6.0, tol);
	EXPECT_NEAR(F.data()[9], 3.0, tol);
	EXPECT_NEAR(F.data()[10], 6.0, tol);
	EXPECT_NEAR(F.data()[11], 9.0, tol);

	// apply natural bcs
	bcApplicator.applyNaturalBCs<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPBdy, HeatEqBundle::QuadratureBoundaryType>(mesh2D, *topoDOF2D, NaturalBCRegistry, t, evalEle, quadBdy, F);

	// tests after applying natural bcs
	EXPECT_NEAR(F.data()[0], 1.0/6.0 - 1.0, tol);
	EXPECT_NEAR(F.data()[1], 1.0/3.0 - 2.0, tol);
	EXPECT_NEAR(F.data()[2], 1.0/2.0 - 3.0, tol);
	EXPECT_NEAR(F.data()[3], 1.0, tol);
	EXPECT_NEAR(F.data()[4], 2.0, tol);
	EXPECT_NEAR(F.data()[5], 3.0, tol);
	EXPECT_NEAR(F.data()[6], 2.0, tol);
	EXPECT_NEAR(F.data()[7], 4.0, tol);
	EXPECT_NEAR(F.data()[8], 6.0, tol);
	EXPECT_NEAR(F.data()[9], 3.0, tol);
	EXPECT_NEAR(F.data()[10], 6.0, tol);
	EXPECT_NEAR(F.data()[11], 9.0, tol);

	// apply essential bcs
	bcApplicator.applyEssentialBCs<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType>(mesh2D, *topoDOF2D, EssentialBCRegistry, t, constantConductivityModel, operatorForms, evalEle, quadVol, F);

	// tests after applying essential bcs
	EXPECT_NEAR(F.data()[0], 1.0/6.0 - 1.0, tol);
	EXPECT_NEAR(F.data()[1], 1.0/3.0 - 2.0, tol);
	EXPECT_NEAR(F.data()[2], 1.0/2.0 - 3.0 + 2.0, tol);
	EXPECT_NEAR(F.data()[3], 1.0, tol);
	EXPECT_NEAR(F.data()[4], 2.0, tol);
	EXPECT_NEAR(F.data()[5], 3.0 + 4.0, tol);
	EXPECT_NEAR(F.data()[6], 2.0, tol);
	EXPECT_NEAR(F.data()[7], 4.0, tol);
	EXPECT_NEAR(F.data()[8], 6.0 + 4.0, tol);
	EXPECT_NEAR(F.data()[9], 3.0 + 1.0, tol);
	EXPECT_NEAR(F.data()[10], 6.0 + 2.0, tol);
	EXPECT_NEAR(F.data()[11], 9.0 + 17.0/3.0, tol);

}
