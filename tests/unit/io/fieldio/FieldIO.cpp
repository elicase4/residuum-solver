#include <filesystem>
#include <fstream>
#include <vector>
#include <gtest/gtest.h>

#include "core/Types.hpp"
#include "equation/heateq/HeatEquation.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/boundary/BoundaryCondition.hpp"
#include "fem/dof/DOFOrdering.hpp"
#include "io/fieldio/FieldIO.hpp"
#include "mesh/ElementFamily.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/generator/BlockMesh2D.hpp"
#include "topology/TopologicalDOF.hpp"

using namespace residuum;

class FieldIOTest : public ::testing::Test {
protected:

	const Real x0 = 0.0;
	const Real x1 = 1.0;
	const Real y0 = 0.0;
	const Real y1 = 1.0;
	const Index nx = 2;
	const Index ny = 2;

	static constexpr Index nsd = 2;
	static constexpr Index Px = 1;
	static constexpr Index Py = 1;
	static constexpr Index numQuadPoint = 2;
	static constexpr Index dofsPerNode = 2;

	using HeatEqBundle = equation::HeatEquation<nsd, 2, mesh::ElementFamily::Quad>;

	HeatEqBundle::Basis basis{Px, Py};

	mesh::generator::BlockMesh2D gen{nx, ny, x0, x1, y0, y1, Px, Py};
	mesh::Mesh mesh;

	std::unique_ptr<topology::TopologicalDOF<dofsPerNode>> topoDOF;
	
	fem::boundary::EssentialBoundaryRegistry EssentialBCRegistry;

	HeatEqBundle::ConductivityModel constantConductivityModel;

	// operator form
	HeatEqBundle::DiffusionForm diffusionForm;
	fem::form::FormRegistry<HeatEqBundle::DiffusionForm> operatorForms{diffusionForm};
	
	static constexpr auto g = [](Real, const Real* x, Real* out){ out[0] = 100.0 + x[0]; out[1] = 200.0 + x[1]; };
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc0;
	
	void SetUp() override {

		mesh = gen.generate();

		constantConductivityModel.setConstant(1.0);
	
		bc0 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{0, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc0);

	}

};

TEST_F(FieldIOTest, reconstructNodalFieldInterleaved){
	
	// dof parameters
	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Interleaved;

	topoDOF = std::make_unique<topology::TopologicalDOF<dofsPerNode>>(mesh, DOFOrdering);
	
	// build constraints
	topoDOF->buildConstraints(basis, EssentialBCRegistry);

	// build free DOF vector
	std::vector<Real> algField(topoDOF->numFreeDOFs());

	for (Index i = 0; i < topoDOF->numFreeDOFs(); ++i) {
		algField[i] = 10.0 + static_cast<Real>(i);
	}

	const auto nodalField = io::fieldio::FieldIO::reconstructNodalField<dofsPerNode>(mesh, *topoDOF, EssentialBCRegistry, 0.0, algField.data());

	ASSERT_EQ(nodalField.size(), mesh.data.numNodes * dofsPerNode);

	for (Index topoIdx = 0; topoIdx < topoDOF->numGlobalDOFs(); ++topoIdx) {

		const Index node = topoDOF->getDOFNode(topoIdx);
		const Index comp = topoIdx - node * dofsPerNode;
		const Real actual = nodalField[node*dofsPerNode + comp];

		if (topoDOF->isConstrained(topoIdx)) {

			const Real* xyz = mesh.getNodeCoord(node);
			const Real expected = (comp == 0) ? (100.0 + xyz[0]) : (200.0 + xyz[1]);
			EXPECT_NEAR(actual, expected, 1e-14);

		} else {

			const Index freeIdx = topoDOF->toAlgebraic(topoIdx);
			EXPECT_NEAR(actual, algField[freeIdx], 1e-14);

		}

	}

}

TEST_F(FieldIOTest, reconstructNodalFieldBlock){

	// dof parameters
	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Block;

	topoDOF = std::make_unique<topology::TopologicalDOF<dofsPerNode>>(mesh, DOFOrdering);
	
	// build constraints
	topoDOF->buildConstraints(basis, EssentialBCRegistry);

	std::vector<Real> algField(topoDOF->numFreeDOFs());

	for (Index i = 0; i < topoDOF->numFreeDOFs(); ++i) {
		algField[i] = 500.0 + static_cast<Real>(i);
	}

	const auto nodalField = io::fieldio::FieldIO::reconstructNodalField<dofsPerNode>(mesh, *topoDOF, EssentialBCRegistry, 0.0, algField.data());

	ASSERT_EQ(nodalField.size(), mesh.data.numNodes * dofsPerNode);

	for (Index topoIdx = 0; topoIdx < topoDOF->numGlobalDOFs(); ++topoIdx) {

		const Index node = topoDOF->getDOFNode(topoIdx);
		const Index comp = topoIdx - node * dofsPerNode;
		const Real actual = nodalField[node*dofsPerNode + comp];

		if (topoDOF->isConstrained(topoIdx)) {

			const Real* xyz = mesh.getNodeCoord(node);
			const Real expected = (comp == 0) ? (100.0 + xyz[0]) : (200.0 + xyz[1]);
			EXPECT_NEAR(actual, expected, 1e-14);

		} else {

			const Index freeIdx = topoDOF->toAlgebraic(topoIdx);
			EXPECT_NEAR(actual, algField[freeIdx], 1e-14);

		}

	}

}

TEST_F(FieldIOTest, WritwVTKContainsFieldNames){

	// dof parameters
	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Block;

	topoDOF = std::make_unique<topology::TopologicalDOF<dofsPerNode>>(mesh, DOFOrdering);
	
	// build constraints
	topoDOF->buildConstraints(basis, EssentialBCRegistry);

	std::vector<Real> algField(topoDOF->numFreeDOFs(), 1.0);

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "field_test.vtk";

	io::fieldio::FieldIO::writeVTK<dofsPerNode>(mesh, *topoDOF, EssentialBCRegistry, 0.0, algField.data(), {"u", "v"}, path.string());

	std::ifstream file(path);
	const std::string content(std::istreambuf_iterator<char>(file), {});

	EXPECT_NE(content.find("POINT_DATA"), std::string::npos);
	EXPECT_NE(content.find("SCALARS u"), std::string::npos);
	EXPECT_NE(content.find("SCALARS v"), std::string::npos);

}

TEST_F(FieldIOTest, writeVTKDOFNameMismatchThrows) {

	// dof parameters
	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Block;

	topoDOF = std::make_unique<topology::TopologicalDOF<dofsPerNode>>(mesh, DOFOrdering);

	// build constraints
	topoDOF->buildConstraints(basis, EssentialBCRegistry);

	std::vector<Real> algField(topoDOF->numFreeDOFs(), 0.0);

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "bad.vtk";

	EXPECT_THROW(io::fieldio::FieldIO::writeVTK<dofsPerNode>(mesh, *topoDOF, EssentialBCRegistry, 0.0, algField.data(), {"u"}, path.string()), std::runtime_error);

}

TEST_F(FieldIOTest, WriteVTUContainsSeparateUnitsAttribute){

	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Block;

	topoDOF = std::make_unique<topology::TopologicalDOF<dofsPerNode>>(mesh, DOFOrdering);
	topoDOF->buildConstraints(basis, EssentialBCRegistry);

	std::vector<Real> algField(topoDOF->numFreeDOFs(), 1.0);

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "field_test.vtu";

	io::fieldio::FieldIO::writeVTU<dofsPerNode>(mesh, *topoDOF, EssentialBCRegistry, 0.0, algField.data(), {"u", "v"}, {"m/s", ""}, path.string());

	std::ifstream file(path);
	const std::string content(std::istreambuf_iterator<char>(file), {});

	EXPECT_NE(content.find("<VTKFile type=\"UnstructuredGrid\""), std::string::npos);
	EXPECT_NE(content.find("Name=\"u\" units=\"m/s\""), std::string::npos);
	// empty unit ("v") must not emit a units attribute at all
	EXPECT_NE(content.find("Name=\"v\" NumberOfComponents"), std::string::npos);
	EXPECT_EQ(content.find("Name=\"v\" units"), std::string::npos);

}

TEST_F(FieldIOTest, writeVTUDOFNameMismatchThrows) {

	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Block;

	topoDOF = std::make_unique<topology::TopologicalDOF<dofsPerNode>>(mesh, DOFOrdering);
	topoDOF->buildConstraints(basis, EssentialBCRegistry);

	std::vector<Real> algField(topoDOF->numFreeDOFs(), 0.0);

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "bad.vtu";

	EXPECT_THROW(io::fieldio::FieldIO::writeVTU<dofsPerNode>(mesh, *topoDOF, EssentialBCRegistry, 0.0, algField.data(), {"u"}, {"", ""}, path.string()), std::runtime_error);

}

TEST_F(FieldIOTest, writeVTUDOFUnitsSizeMismatchThrows) {

	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Block;

	topoDOF = std::make_unique<topology::TopologicalDOF<dofsPerNode>>(mesh, DOFOrdering);
	topoDOF->buildConstraints(basis, EssentialBCRegistry);

	std::vector<Real> algField(topoDOF->numFreeDOFs(), 0.0);

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "bad_units.vtu";

	// dofUnits has 1 entry, dofNames has 2 -- must be rejected
	EXPECT_THROW(io::fieldio::FieldIO::writeVTU<dofsPerNode>(mesh, *topoDOF, EssentialBCRegistry, 0.0, algField.data(), {"u", "v"}, {"m"}, path.string()), std::runtime_error);

}

TEST_F(FieldIOTest, BinaryWriteReadRoundTrip){

	// dof parameters
	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Interleaved;

	topoDOF = std::make_unique<topology::TopologicalDOF<dofsPerNode>>(mesh, DOFOrdering);

	// build constraints
	topoDOF->buildConstraints(basis, EssentialBCRegistry);

	std::vector<Real> algField(topoDOF->numFreeDOFs());

	for (Index i = 0; i < topoDOF->numFreeDOFs(); ++i) {
		algField[i] = 10.0 + static_cast<Real>(i);
	}

	const auto expected = io::fieldio::FieldIO::reconstructNodalField<dofsPerNode>(mesh, *topoDOF, EssentialBCRegistry, 0.0, algField.data());

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "field_test.pndf";

	io::fieldio::FieldIO::writeBinary<dofsPerNode>(mesh, *topoDOF, EssentialBCRegistry, 0.0, algField.data(), path.string());

	const auto loaded = io::fieldio::FieldIO::readBinary<dofsPerNode>(mesh, path.string());

	EXPECT_NEAR(loaded.time, 0.0, 1e-14);

	ASSERT_EQ(loaded.values.size(), expected.size());

	for (Index i = 0; i < expected.size(); ++i){
		EXPECT_NEAR(loaded.values[i], expected[i], 1e-12);
	}

}

TEST_F(FieldIOTest, BinaryWriteRawReadRoundTrip){

	std::vector<Real> nodalField(mesh.data.numNodes * dofsPerNode);
	for (Index i = 0; i < nodalField.size(); ++i) {
		nodalField[i] = 42.0 + static_cast<Real>(i);
	}

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "field_test_raw.pndf";

	io::fieldio::FieldIO::writeBinaryRaw<dofsPerNode>(mesh, 3.5, nodalField, path.string());

	const auto loaded = io::fieldio::FieldIO::readBinary<dofsPerNode>(mesh, path.string());

	EXPECT_NEAR(loaded.time, 3.5, 1e-14);
	ASSERT_EQ(loaded.values.size(), nodalField.size());

	for (Index i = 0; i < nodalField.size(); ++i){
		EXPECT_NEAR(loaded.values[i], nodalField[i], 1e-12);
	}

}

TEST_F(FieldIOTest, BinaryWriteRawSizeMismatchThrows){

	std::vector<Real> wrongSized(mesh.data.numNodes * dofsPerNode - 1);

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "field_test_raw_badsize.pndf";

	EXPECT_THROW(io::fieldio::FieldIO::writeBinaryRaw<dofsPerNode>(mesh, 0.0, wrongSized, path.string()), std::runtime_error);

}

TEST_F(FieldIOTest, BinaryReadNumNodesMismatchThrows){

	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Interleaved;

	topoDOF = std::make_unique<topology::TopologicalDOF<dofsPerNode>>(mesh, DOFOrdering);
	topoDOF->buildConstraints(basis, EssentialBCRegistry);

	std::vector<Real> algField(topoDOF->numFreeDOFs(), 1.0);

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "field_test_numnodesmismatch.pndf";

	io::fieldio::FieldIO::writeBinary<dofsPerNode>(mesh, *topoDOF, EssentialBCRegistry, 0.0, algField.data(), path.string());

	// a differently-refined mesh has a different numNodes -- must be rejected before ever touching meshTag
	mesh::generator::BlockMesh2D genFiner{nx + 1, ny + 1, x0, x1, y0, y1, Px, Py};
	mesh::Mesh finerMesh = genFiner.generate();

	EXPECT_THROW(io::fieldio::FieldIO::readBinary<dofsPerNode>(finerMesh, path.string()), std::runtime_error);

}

TEST_F(FieldIOTest, BinaryReadMeshTagMismatchThrows){

	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Interleaved;

	topoDOF = std::make_unique<topology::TopologicalDOF<dofsPerNode>>(mesh, DOFOrdering);
	topoDOF->buildConstraints(basis, EssentialBCRegistry);

	std::vector<Real> algField(topoDOF->numFreeDOFs(), 1.0);

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "field_test_tagmismatch.pndf";

	io::fieldio::FieldIO::writeBinary<dofsPerNode>(mesh, *topoDOF, EssentialBCRegistry, 0.0, algField.data(), path.string());

	// same numNodes/numElements as `mesh` (identical nx/ny/Px/Py), different geometry -- the numNodes check alone can't catch this, only computeMeshTag can
	mesh::generator::BlockMesh2D genOther{nx, ny, x0, x1 + 1.0, y0, y1, Px, Py};
	mesh::Mesh otherMesh = genOther.generate();

	ASSERT_EQ(otherMesh.data.numNodes, mesh.data.numNodes);

	EXPECT_THROW(io::fieldio::FieldIO::readBinary<dofsPerNode>(otherMesh, path.string()), std::runtime_error);

}
