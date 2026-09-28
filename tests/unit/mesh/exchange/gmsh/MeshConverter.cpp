#include <filesystem>
#include <gtest/gtest.h>

#include "core/Types.hpp"
#include "io/GmshReader.hpp"
#include "io/MeshIO.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/exchange/gmsh/IntermediateMesh.hpp"
#include "mesh/exchange/gmsh/MeshConverter.hpp"

using namespace residuum;

TEST(MeshConverter, QuadP1Structured){

	mesh::exchange::gmsh::IntermediateMesh gmshMesh;

	const std::filesystem::path input_mesh_path = std::filesystem::path(TEST_DATA_PATH) / "mesh/gmsh/msh/quad_p1.msh";

	io::GmshReader::read(gmshMesh, input_mesh_path.string());

	mesh::Mesh mesh;

	std::unordered_map<Int, Int> boundaryMap = {{1,2}, {2,1}, {3,3}, {4,0}};

	mesh::exchange::gmsh::MeshConverter::toSolverMesh(mesh, gmshMesh, boundaryMap);

	// check metadata
	EXPECT_TRUE(mesh.isValid());
	EXPECT_EQ(mesh.data.parametricDim, 2);
	EXPECT_EQ(mesh.data.spatialDim, 2);
	EXPECT_EQ(mesh.data.numNodes, 9);
	EXPECT_EQ(mesh.data.numElements, 4);
	EXPECT_EQ(mesh.data.nodesPerElement, 4);
	EXPECT_EQ(mesh.data.facesPerElement, 4);
	EXPECT_EQ(mesh.data.basisOrder, std::vector<Index>({1,1}));
	
	// check coordinates
	std::set<std::pair<Real,Real>> expected = {
		{0.0, 0.0},
		{0.5, 0.0},
		{1.0, 0.0},
		{0.0, 0.5},
		{0.5, 0.5},
		{1.0, 0.5},
		{0.0, 1.0},
		{0.5, 1.0},
		{1.0, 1.0}
	};

	std::set<std::pair<Real,Real>> found;
	for (Index n = 0; n < mesh.data.numNodes; ++n){
		found.insert({mesh.data.xyz[n*2+0], mesh.data.xyz[n*2+1]});
	}
	EXPECT_EQ(found, expected);

	// check connectivity (local canonical ordering)
	for (Index e = 0; e < mesh.data.numElements; ++e) {

		const Index* ien = mesh.getElementNodes(e);
		
		auto xy = [&](Index n){
			return std::pair<Real,Real>{mesh.data.xyz[2*n+0], mesh.data.xyz[2*n+1]};
		};

		auto p0 = xy(ien[0]);
		auto p1 = xy(ien[1]);
		auto p2 = xy(ien[2]);
		auto p3 = xy(ien[3]);

		// xi ordering
		EXPECT_LT(p0.first, p1.first);
		EXPECT_LT(p2.first, p3.first);

		// eta ordering
		EXPECT_LT(p0.second, p2.second);
		EXPECT_LT(p1.second, p3.second);

	}

	// check boundary tags (local canonical ordering)
	const Real tol = 1e-12;
	for (Index e = 0; e < mesh.data.numElements; ++e) {
		
		const Index* ien = mesh.getElementNodes(e);
		Real xc = 0.0;
		Real yc = 0.0;

		for (Index a = 0; a < mesh.data.nodesPerElement; ++a){
			xc += mesh.data.xyz[2*ien[a]+0];
			yc += mesh.data.xyz[2*ien[a]+1];
		}

		xc /= mesh.data.nodesPerElement;
		yc /= mesh.data.nodesPerElement;

		Int* rng = mesh.getBoundaryTag(e);

		Int expectedRNG[4] = {
			(xc < tol + 0.25) ? 0 : -1,
			(xc > 0.75 - tol) ? 1 : -1,
			(yc < tol + 0.25) ? 2 : -1,
			(yc > 0.75 - tol) ? 3 : -1,
		};

		for (Index f = 0; f < 4; ++f){
			EXPECT_EQ(rng[f], expectedRNG[f]);
		}

	}

	// write mesh vtk file
	const std::filesystem::path output_path = std::filesystem::path(TEST_OUTPUT_PATH) / "quad_p1_mesh.vtk";
	residuum::io::MeshIO::writeVTK(mesh, output_path.string());

}

TEST(MeshConverter, HexP1Structured){
	
	mesh::exchange::gmsh::IntermediateMesh gmshMesh;

	const std::filesystem::path input_mesh_path = std::filesystem::path(TEST_DATA_PATH) / "mesh/gmsh/msh/hex_p1.msh";

	io::GmshReader::read(gmshMesh, input_mesh_path.string());

	mesh::Mesh mesh;

	std::unordered_map<Int, Int> boundaryMap = {{1,0}, {2,1}, {3,2}, {4,3}, {5,4}, {6,5}};

	mesh::exchange::gmsh::MeshConverter::toSolverMesh(mesh, gmshMesh, boundaryMap);

	// check metadata
	EXPECT_TRUE(mesh.isValid());
	EXPECT_EQ(mesh.data.parametricDim, 3);
	EXPECT_EQ(mesh.data.spatialDim, 3);
	EXPECT_EQ(mesh.data.numNodes, 27);
	EXPECT_EQ(mesh.data.numElements, 8);
	EXPECT_EQ(mesh.data.nodesPerElement, 8);
	EXPECT_EQ(mesh.data.facesPerElement, 6);
	EXPECT_EQ(mesh.data.basisOrder, std::vector<Index>({1,1,1}));
	
	// check coordinates
	std::set<std::tuple<Real,Real,Real>> expected;
	for (Index k = 0; k <= 2; ++k){
		for (Index j = 0; j <= 2; ++j){
			for (Index i = 0; i <= 2; ++i){
				expected.insert({0.5*i, 0.5*j, 0.5*k});
			}
		}
	}

	std::set<std::tuple<Real,Real,Real>> found;
	for (Index n = 0; n < mesh.data.numNodes; ++n){
		found.insert({mesh.data.xyz[n*3+0], mesh.data.xyz[n*3+1], mesh.data.xyz[3*n+2]});
	}
	EXPECT_EQ(found, expected);

	// check connectivity (local canonical ordering)
	for (Index e = 0; e < mesh.data.numElements; ++e) {

		const Index* ien = mesh.getElementNodes(e);
		
		auto xy = [&](Index n){
			return std::tuple<Real,Real,Real>{mesh.data.xyz[3*n+0], mesh.data.xyz[3*n+1], mesh.data.xyz[3*n+2]};
		};

		auto p0 = xy(ien[0]);
		auto p1 = xy(ien[1]);
		auto p2 = xy(ien[2]);
		auto p3 = xy(ien[3]);
		auto p4 = xy(ien[4]);
		auto p5 = xy(ien[5]);
		auto p6 = xy(ien[6]);
		auto p7 = xy(ien[7]);

		// xi ordering
		EXPECT_LT(std::get<0>(p0), std::get<0>(p1));
		EXPECT_LT(std::get<0>(p2), std::get<0>(p3));
		EXPECT_LT(std::get<0>(p4), std::get<0>(p5));
		EXPECT_LT(std::get<0>(p6), std::get<0>(p7));

		// eta ordering
		EXPECT_LT(std::get<1>(p0), std::get<1>(p2));
		EXPECT_LT(std::get<1>(p1), std::get<1>(p3));
		EXPECT_LT(std::get<1>(p4), std::get<1>(p6));
		EXPECT_LT(std::get<1>(p5), std::get<1>(p7));

		// zeta ordering
		EXPECT_LT(std::get<2>(p0), std::get<2>(p4));
		EXPECT_LT(std::get<2>(p1), std::get<2>(p5));
		EXPECT_LT(std::get<2>(p2), std::get<2>(p6));
		EXPECT_LT(std::get<2>(p3), std::get<2>(p7));

	}

	// check boundary tags (local canonical ordering)
	const Real tol = 1e-12;
	for (Index e = 0; e < mesh.data.numElements; ++e) {
		
		const Index* ien = mesh.getElementNodes(e);
		Real xc = 0.0;
		Real yc = 0.0;
		Real zc = 0.0;

		for (Index a = 0; a < mesh.data.nodesPerElement; ++a){
			xc += mesh.data.xyz[3*ien[a]+0];
			yc += mesh.data.xyz[3*ien[a]+1];
			zc += mesh.data.xyz[3*ien[a]+2];
		}

		xc /= mesh.data.nodesPerElement;
		yc /= mesh.data.nodesPerElement;
		zc /= mesh.data.nodesPerElement;

		Int* rng = mesh.getBoundaryTag(e);

		Int expectedRNG[6] = {
			(xc < tol + 0.25) ? 0 : -1,
			(xc > 0.75 - tol) ? 1 : -1,
			(yc < tol + 0.25) ? 2 : -1,
			(yc > 0.75 - tol) ? 3 : -1,
			(zc < tol + 0.25) ? 4 : -1,
			(zc > 0.75 - tol) ? 5 : -1
		};

		for (Index f = 0; f < 6; ++f){
			EXPECT_EQ(rng[f], expectedRNG[f]);
		}

	}

	// write mesh vtk file
	const std::filesystem::path output_path = std::filesystem::path(TEST_OUTPUT_PATH) / "hex_p1_mesh.vtk";
	residuum::io::MeshIO::writeVTK(mesh, output_path.string());

}
