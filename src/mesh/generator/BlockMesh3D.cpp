#include "mesh/generator/BlockMesh3D.hpp"

#include <stdexcept>
#include <unordered_map>

residuum::mesh::Mesh residuum::mesh::generator::BlockMesh3D::generate(const std::unordered_map<Int, Int>&) const {

	if (nx_ <= 0 || ny_ <= 0 || nz_ <= 0 || px_ <= 0 || py_ <= 0 || pz_ <= 0) {
		throw std::invalid_argument("BlockMesh3D: nx, ny, nz, px, py, pz must all be positive");
	}

	Mesh mesh;

	initializeData(mesh.data);
	generateNodes(mesh.data);
	generateElements(mesh.data);
	generateBoundaryTags(mesh.data);

	mesh.data.elementFamily = ElementFamily::Hex;
	mesh.data.basisType = BasisType::Lagrange;

	return mesh;

}

void residuum::mesh::generator::BlockMesh3D::initializeData(Data& data) const {

	// dimensions
	data.parametricDim = 3;
	data.spatialDim = 3;

	// basis function info
	data.basisOrder = {px_, py_, pz_};
	data.nodesPerElement = (px_ + 1)*(py_ + 1)*(pz_ + 1);

	// mesh size info
	data.numNodes = (nx_*px_ + 1) * (ny_*py_ + 1) * (nz_*pz_ + 1);
	data.numElements = nx_ * ny_ * nz_;

	// boundary info
	data.facesPerElement = 6;

	// size arrays
	data.xyz.resize(data.numNodes * data.spatialDim);
	data.ien.resize(data.numElements * data.nodesPerElement);
	data.rng.resize(data.numElements * data.facesPerElement);

}

void residuum::mesh::generator::BlockMesh3D::generateNodes(Data& data) const {

	const Index NX = nx_ * px_ + 1;
	const Index NY = ny_ * py_ + 1;
	const Index NZ = nz_ * pz_ + 1;

	const Real dx = (x1_ - x0_) / static_cast<Real>(nx_ * px_);
	const Real dy = (y1_ - y0_) / static_cast<Real>(ny_ * py_);
	const Real dz = (z1_ - z0_) / static_cast<Real>(nz_ * pz_);

	Index node = 0;

	for (Index k = 0; k < NZ; ++k) {
		for (Index j = 0; j < NY; ++j) {
			for (Index i = 0; i < NX; ++i) {

				data.xyz[data.spatialDim * node] = x0_ + static_cast<Real>(i) * dx;
				data.xyz[data.spatialDim * node + 1] = y0_ + static_cast<Real>(j) * dy;
				data.xyz[data.spatialDim * node + 2] = z0_ + static_cast<Real>(k) * dz;

				node++;

			}
		}
	}

}

void residuum::mesh::generator::BlockMesh3D::generateElements(Data& data) const {

	const Index NX = nx_ * px_ + 1;
	const Index NY = ny_ * py_ + 1;

	// populate ien array
	for (Index ele_z = 0; ele_z < nz_; ++ele_z){
		for (Index ele_y = 0; ele_y < ny_; ++ele_y){
			for (Index ele_x = 0; ele_x < nx_; ++ele_x){

				Index ele = ele_z * (nx_ * ny_) + ele_y * nx_ + ele_x;

				for (Index a = 0; a < data.nodesPerElement; ++a){

					Index a_x = a % (px_ + 1);
					Index a_y = (a / (px_ + 1)) % (py_ + 1);
					Index a_z = a / ((px_ + 1)*(py_ + 1));

					Index node_x = ele_x * px_ + a_x;
					Index node_y = ele_y * py_ + a_y;
					Index node_z = ele_z * pz_ + a_z;

					data.ien[ele * data.nodesPerElement + a] = node_z * (NY * NX) + node_y * NX + node_x;

				}

			}
		}
	}

}

void residuum::mesh::generator::BlockMesh3D::generateBoundaryTags(Data& data) const {

	// boundary tags
	static constexpr Int LEFT   = 0;
	static constexpr Int RIGHT  = 1;
	static constexpr Int FRONT  = 2;
	static constexpr Int BACK   = 3;
	static constexpr Int BOTTOM = 4;
	static constexpr Int TOP    = 5;

	// populate rng array
	for (Index ele_z = 0; ele_z < nz_; ++ele_z){
		for (Index ele_y = 0; ele_y < ny_; ++ele_y){
			for (Index ele_x = 0; ele_x < nx_; ++ele_x){

				// get flat id
				Index ele = ele_z * (nx_ * ny_) + ele_y * nx_ + ele_x;

				// set boundary tags
				data.rng[ele * data.facesPerElement + LEFT]   = (ele_x == 0)        ? LEFT   : -1;
				data.rng[ele * data.facesPerElement + RIGHT]  = (ele_x == (nx_ - 1)) ? RIGHT  : -1;
				data.rng[ele * data.facesPerElement + FRONT]  = (ele_y == 0)        ? FRONT  : -1;
				data.rng[ele * data.facesPerElement + BACK]   = (ele_y == (ny_ - 1)) ? BACK   : -1;
				data.rng[ele * data.facesPerElement + BOTTOM] = (ele_z == 0)        ? BOTTOM : -1;
				data.rng[ele * data.facesPerElement + TOP]    = (ele_z == (nz_ - 1)) ? TOP    : -1;

			}
		}
	}
}
