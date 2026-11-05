#ifndef RESIDUUM_MESH_GENERATOR_BLOCKMESH3D_HPP
#define RESIDUUM_MESH_GENERATOR_BLOCKMESH3D_HPP

#include "mesh/generator/MeshGenerator.hpp"

namespace residuum {

	namespace mesh {

		namespace generator {

			class BlockMesh3D : public MeshGenerator {
			public:
				BlockMesh3D(Index nx_, Index ny_, Index nz_, Real x0_, Real x1_, Real y0_, Real y1_, Real z0_, Real z1_, Index px_, Index py_, Index pz_): nx_(nx_), ny_(ny_), nz_(nz_), x0_(x0_), x1_(x1_), y0_(y0_), y1_(y1_), z0_(z0_), z1_(z1_), px_(px_), py_(py_), pz_(pz_) {};

				Mesh generate(const std::unordered_map<Int, Int>& physicalGroupMap = {}) const override;

			private:

				void initializeData(Data& data) const;

				void generateNodes(Data& data) const;

				void generateElements(Data& data) const;

				void generateBoundaryTags(Data& data) const;

				Index nx_;
				Index ny_;
				Index nz_;

				Real x0_;
				Real x1_;
				Real y0_;
				Real y1_;
				Real z0_;
				Real z1_;

				Index px_;
				Index py_;
				Index pz_;

			}; // class BlockMesh3D

		} // namespace generator

	} // namespace mesh

} // namespace residuum

#endif
