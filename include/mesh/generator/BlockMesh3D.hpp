#ifndef RESIDUUM_MESH_GENERATOR_BLOCKMESH3D_HPP
#define RESIDUUM_MESH_GENERATOR_BLOCKMESH3D_HPP

#include "mesh/generator/MeshGenerator.hpp"

namespace residuum {

	namespace mesh {

		namespace generator {

			class BlockMesh3D : public MeshGenerator {
			public:
				BlockMesh3D(Index nx_, Index ny_, Index nz_, Real x0_, Real x1_, Real y0_, Real y1_, Real z0_, Real z1_, Index px_, Index py_, Index pz_): nx(nx_), ny(ny_), nz(nz_), x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), px(px_), py(py_), pz(pz_) {};

				Mesh generate(const std::unordered_map<Int, Int>& physicalGroupMap = {}) const override;

			private:

				void initializeData(Data& data) const;

				void generateNodes(Data& data) const;

				void generateElements(Data& data) const;

				void generateBoundaryTags(Data& data) const;

				Index nx;
				Index ny;
				Index nz;

				Real x0;
				Real x1;
				Real y0;
				Real y1;
				Real z0;
				Real z1;

				Index px;
				Index py;
				Index pz;

			}; // class BlockMesh3D

		} // namespace generator

	} // namespace mesh

} // namespace residuum

#endif
