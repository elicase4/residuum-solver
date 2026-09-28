#include <filesystem>
#include <gtest/gtest.h>

#include "core/FEM.hpp"
#include "core/IO.hpp"
#include "core/LinAlg.hpp"
#include "core/Mesh.hpp"
#include "core/Topology.hpp"
#include "core/Types.hpp"

#include "equation/heateq/HeatEquation.hpp"
#include "mesh/ElementFamily.hpp"
#include "mesh/generator/BlockMesh2D.hpp"

using namespace residuum;

namespace {

	constexpr auto sourceExpr = [](Real, const Real* x, Real* out) { out[0] = x[0] - x[1]; };
	constexpr auto fluxExpr = [](Real, const Real* x, Real* out) { out[0] = x[0]; out[1] = x[1]; };

} // namespace

class NodalDataFormEquivalence : public ::testing::Test {
protected:

	static constexpr Index nsd = 2;
	static constexpr Index Px = 1;
	static constexpr Index Py = 1;
	static constexpr Index numQuadPoint = 2;

	using BackendType = linalg::types::backend::CPU;
	using HeatEqBundle = equation::HeatEquation<nsd, 2, mesh::ElementFamily::Quad>;

	HeatEqBundle::Basis basis{Px, Py};
	HeatEqBundle::QuadratureVolumeType quadVol{numQuadPoint, numQuadPoint};
	HeatEqBundle::QuadratureBoundaryType quadBdy{numQuadPoint};
	HeatEqBundle::EvalEle evalEle{basis};

	mesh::generator::BlockMesh2D gen{4, 4, -1.0, 1.0, -1.0, 1.0, Px, Py};
	mesh::Mesh mesh;
	std::unique_ptr<topology::TopologicalDOF<HeatEqBundle::NumDOFs>> topoDOF;

	fem::boundary::EssentialBoundaryRegistry essentialBCs;

	fem::assembly::Assembler<BackendType> assembler;
	fem::boundary::BoundaryApplicator<BackendType> bcApplicator;

	HeatEqBundle::DefaultModel defaultModel;
	HeatEqBundle::DefaultModelBdy defaultModelBdy;

	void SetUp() override {

		mesh = gen.generate();
		topoDOF = std::make_unique<topology::TopologicalDOF<HeatEqBundle::NumDOFs>>(mesh, fem::dof::DOFOrdering::Interleaved);

		topoDOF->buildConstraints(basis, essentialBCs);

	}

};

TEST_F(NodalDataFormEquivalence, NodalSourceFormMatchesExpressionSourceFormForLinearField) {

	using ExpressionSourceFormT = HeatEqBundle::SourceForm<decltype(sourceExpr)>;
	fem::form::FormRegistry<ExpressionSourceFormT> exprForms{sourceExpr};

	// write the same field the expression evaluates, node by node, as a .pndf fixture
	std::vector<Real> sourceField(mesh.data.numNodes * HeatEqBundle::NumDOFs);
	for (Index nodeID = 0; nodeID < mesh.data.numNodes; ++nodeID) {
		sourceExpr(0.0, mesh.getNodeCoord(nodeID), &sourceField[nodeID * HeatEqBundle::NumDOFs]);
	}
	const auto sourcePath = (std::filesystem::path(TEST_OUTPUT_PATH) / "source_field.pndf").string();
	io::fieldio::FieldIO::writeBinaryRaw<HeatEqBundle::NumDOFs>(mesh, 0.0, sourceField, sourcePath);

	using NodalSourceSourceT = io::fieldio::NodalFileValueSource<HeatEqBundle::NumDOFs>;
	using NodalSourceFormT = HeatEqBundle::NodalSourceForm<NodalSourceSourceT>;
	NodalSourceSourceT nodalSource(mesh, sourcePath);
	fem::form::FormRegistry<NodalSourceFormT> nodalForms{nodalSource};

	auto U = assembler.createVector<HeatEqBundle::NumDOFs>(mesh, *topoDOF);
	auto Fexpr = assembler.createVector<HeatEqBundle::NumDOFs>(mesh, *topoDOF);
	auto Ffile = assembler.createVector<HeatEqBundle::NumDOFs>(mesh, *topoDOF);
	U.zero();
	Fexpr.zero();
	Ffile.zero();

	assembler.assembleVector<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::DefaultModel, decltype(exprForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh, *topoDOF, 0.0, defaultModel, exprForms, evalEle, quadVol, U, nullptr, {nullptr}, Fexpr, nullptr);
	assembler.assembleVector<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::DefaultModel, decltype(nodalForms), HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh, *topoDOF, 0.0, defaultModel, nodalForms, evalEle, quadVol, U, nullptr, {nullptr}, Ffile, nullptr);

	ASSERT_EQ(Fexpr.size(), Ffile.size());
	for (Index i = 0; i < Fexpr.size(); ++i) {
		EXPECT_NEAR(Ffile.data()[i], Fexpr.data()[i], 1e-10);
	}

}

TEST_F(NodalDataFormEquivalence, NodalFluxFormMatchesExpressionFluxFormForLinearField) {

	static constexpr Int boundaryTag = 0;

	using FluxFunctionExpressionT = HeatEqBundle::FluxBC<decltype(fluxExpr)>;
	using ExpressionFluxFormT = HeatEqBundle::FluxForm<decltype(fluxExpr)>;
	fem::form::FormRegistry<ExpressionFluxFormT> exprForms{fluxExpr};

	auto bcExpr = std::shared_ptr<fem::boundary::BoundaryCondition<FluxFunctionExpressionT>>(new fem::boundary::BoundaryCondition<FluxFunctionExpressionT>{boundaryTag, {fem::boundary::BCCategory::Natural}, FluxFunctionExpressionT{fluxExpr}});
	fem::boundary::NaturalBoundaryRegistry<HeatEqBundle::EvalQPBdy> exprRegistry;
	exprRegistry.registerBC<FluxFunctionExpressionT>(bcExpr, exprForms, defaultModelBdy);

	// write the same vector field the expression evaluates, node by node, as a .pndf fixture
	std::vector<Real> fluxField(mesh.data.numNodes * HeatEqBundle::NumDOFs * HeatEqBundle::SpatialDim);
	for (Index nodeID = 0; nodeID < mesh.data.numNodes; ++nodeID) {
		fluxExpr(0.0, mesh.getNodeCoord(nodeID), &fluxField[nodeID * HeatEqBundle::NumDOFs * HeatEqBundle::SpatialDim]);
	}
	const auto fluxPath = (std::filesystem::path(TEST_OUTPUT_PATH) / "flux_field.pndf").string();
	io::fieldio::FieldIO::writeBinaryRaw<HeatEqBundle::NumDOFs * HeatEqBundle::SpatialDim>(mesh, 0.0, fluxField, fluxPath);

	using NodalFluxSourceT = io::fieldio::NodalFileValueSource<HeatEqBundle::NumDOFs * HeatEqBundle::SpatialDim>;
	using NodalFluxFormT = HeatEqBundle::NodalFluxForm<NodalFluxSourceT>;
	using FluxFunctionFileT = HeatEqBundle::FluxBC<io::fieldio::NodalValueSourceAdapter<NodalFluxSourceT>>;

	NodalFluxSourceT nodalFluxSource(mesh, fluxPath);
	fem::form::FormRegistry<NodalFluxFormT> nodalForms{nodalFluxSource};

	auto bcFile = std::shared_ptr<fem::boundary::BoundaryCondition<FluxFunctionFileT>>(new fem::boundary::BoundaryCondition<FluxFunctionFileT>{boundaryTag, {fem::boundary::BCCategory::Natural}, FluxFunctionFileT{mesh, fluxPath}});
	fem::boundary::NaturalBoundaryRegistry<HeatEqBundle::EvalQPBdy> fileRegistry;
	fileRegistry.registerBC<FluxFunctionFileT>(bcFile, nodalForms, defaultModelBdy);

	auto Fexpr = assembler.createVector<HeatEqBundle::NumDOFs>(mesh, *topoDOF);
	auto Ffile = assembler.createVector<HeatEqBundle::NumDOFs>(mesh, *topoDOF);
	Fexpr.zero();
	Ffile.zero();

	bcApplicator.applyNaturalBCs<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPBdy, HeatEqBundle::QuadratureBoundaryType>(mesh, *topoDOF, exprRegistry, 0.0, evalEle, quadBdy, Fexpr);
	bcApplicator.applyNaturalBCs<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPBdy, HeatEqBundle::QuadratureBoundaryType>(mesh, *topoDOF, fileRegistry, 0.0, evalEle, quadBdy, Ffile);

	ASSERT_EQ(Fexpr.size(), Ffile.size());
	for (Index i = 0; i < Fexpr.size(); ++i) {
		EXPECT_NEAR(Ffile.data()[i], Fexpr.data()[i], 1e-10);
	}

}
