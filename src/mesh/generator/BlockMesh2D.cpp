#include "mesh/generator/BlockMesh2D.hpp"

#include <stdexcept>
#include <unordered_map>

residuum::mesh::Mesh residuum::mesh::generator::BlockMesh2D::generate(const std::unordered_map<Int, Int>&) const {

	if (nx_ <= 0 || ny_ <= 0 || px_ <= 0 || py_ <= 0) {
		throw std::invalid_argument("BlockMesh2D: nx, ny, px, py must all be positive");
	}

	Mesh mesh;

	initializeData(mesh.data);
	generateNodes(mesh.data);
	generateElements(mesh.data);
	generateBoundaryTags(mesh.data);

	mesh.data.elementFamily = ElementFamily::Quad;
	mesh.data.basisType = BasisType::Lagrange;

	return mesh;

}

void residuum::mesh::generator::BlockMesh2D::initializeData(Data& data) const {

	// dimensions
	data.parametricDim = 2;
	data.spatialDim = 2;

	// basis function info
	data.basisOrder = {px_, py_};
	data.nodesPerElement = (px_ + 1)*(py_ + 1);

	// mesh size info
	data.numNodes = (nx_*px_ + 1) * (ny_*py_ + 1);
	data.numElements = nx_ * ny_;

	// boundary info
	data.facesPerElement = 4;

	// size arrays
	data.xyz.resize(data.numNodes * data.spatialDim);
	data.ien.resize(data.numElements * data.nodesPerElement);
	data.rng.resize(data.numElements * data.facesPerElement);

}

void residuum::mesh::generator::BlockMesh2D::generateNodes(Data& data) const {

	const Index NX = nx_ * px_ + 1;
	const Index NY = ny_ * py_ + 1;

	const Real dx = (x1_ - x0_) / static_cast<Real>(nx_ * px_);
	const Real dy = (y1_ - y0_) / static_cast<Real>(ny_ * py_);

	Index node = 0;

	for (Index j = 0; j < NY; ++j) {
		for (Index i = 0; i < NX; ++i) {

			data.xyz[data.spatialDim * node] = x0_ + static_cast<Real>(i) * dx;
			data.xyz[data.spatialDim * node + 1] = y0_ + static_cast<Real>(j) * dy;

			node++;

		}
	}

}

void residuum::mesh::generator::BlockMesh2D::generateElements(Data& data) const {

	// populate ien array
	for (Index ele_y = 0; ele_y < ny_; ++ele_y){
		for (Index ele_x = 0; ele_x < nx_; ++ele_x){

			Index ele = ele_y * nx_ + ele_x;

			for (Index a = 0; a < data.nodesPerElement; ++a){

				Index a_x  = a % (px_ + 1);
				Index a_y  = a / (px_ + 1);

				Index node_x = ele_x * px_ + a_x;
				Index node_y = ele_y * py_ + a_y;

				data.ien[ele * data.nodesPerElement + a] = node_y * (nx_ * px_ + 1) + node_x;

			}

		}
	}

}

void residuum::mesh::generator::BlockMesh2D::generateBoundaryTags(Data& data) const {

	// boundary tags
	static constexpr Int LEFT   = 0;
	static constexpr Int RIGHT  = 1;
	static constexpr Int BOTTOM = 2;
	static constexpr Int TOP    = 3;

	// populate rng array
	for (Index ele_y = 0; ele_y < ny_; ++ele_y){
		for (Index ele_x = 0; ele_x < nx_; ++ele_x){

			// get flat id
			Index ele = ele_y * nx_ + ele_x;

			// set boundary tags
			data.rng[ele * data.facesPerElement + LEFT]   = (ele_x == 0)        ? LEFT   : -1;
			data.rng[ele * data.facesPerElement + RIGHT]  = (ele_x == (nx_ - 1)) ? RIGHT  : -1;
			data.rng[ele * data.facesPerElement + BOTTOM] = (ele_y == 0)        ? BOTTOM : -1;
			data.rng[ele * data.facesPerElement + TOP]    = (ele_y == (ny_ - 1)) ? TOP    : -1;

		}
	}
}
