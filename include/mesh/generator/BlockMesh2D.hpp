#ifndef RESIDUUM_MESH_GENERATOR_BLOCKMESH2D_HPP
#define RESIDUUM_MESH_GENERATOR_BLOCKMESH2D_HPP

#include "mesh/generator/MeshGenerator.hpp"

namespace residuum {

	namespace mesh {

		namespace generator {

			class BlockMesh2D : public MeshGenerator {
			public:
				BlockMesh2D(Index nx_, Index ny_, Real x0_, Real x1_, Real y0_, Real y1_, Index px_, Index py_): nx(nx_), ny(ny_), x0(x0_), x1(x1_), y0(y0_), y1(y1_), px(px_), py(py_) {};

				Mesh generate(const std::unordered_map<Int, Int>& physicalGroupMap = {}) const override;

			private:

				void initializeData(Data& data) const;

				void generateNodes(Data& data) const;

				void generateElements(Data& data) const;

				void generateBoundaryTags(Data& data) const;

				Index nx;
				Index ny;

				Real x0;
				Real x1;
				Real y0;
				Real y1;

				Index px;
				Index py;

			}; // class BlockMesh2D

		} // namespace generator

	} // namespace mesh

} // namespace residuum

#endif
