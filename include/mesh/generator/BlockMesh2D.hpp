#ifndef RESIDUUM_MESH_GENERATOR_BLOCKMESH2D_HPP
#define RESIDUUM_MESH_GENERATOR_BLOCKMESH2D_HPP

#include "mesh/generator/MeshGenerator.hpp"

namespace residuum {

	namespace mesh {

		namespace generator {

			class BlockMesh2D : public MeshGenerator {
			public:
				BlockMesh2D(Index nx_, Index ny_, Real x0_, Real x1_, Real y0_, Real y1_, Index px_, Index py_): nx_(nx_), ny_(ny_), x0_(x0_), x1_(x1_), y0_(y0_), y1_(y1_), px_(px_), py_(py_) {};

				Mesh generate(const std::unordered_map<Int, Int>& physicalGroupMap = {}) const override;

			private:

				void initializeData(Data& data) const;

				void generateNodes(Data& data) const;

				void generateElements(Data& data) const;

				void generateBoundaryTags(Data& data) const;

				Index nx_;
				Index ny_;

				Real x0_;
				Real x1_;
				Real y0_;
				Real y1_;

				Index px_;
				Index py_;

			}; // class BlockMesh2D

		} // namespace generator

	} // namespace mesh

} // namespace residuum

#endif
