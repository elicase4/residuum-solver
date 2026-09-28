#ifndef RESIDUUM_MESH_DATA_HPP
#define RESIDUUM_MESH_DATA_HPP

#include "core/Types.hpp"
#include "mesh/BasisType.hpp"
#include "mesh/ElementFamily.hpp"
#include <vector>

namespace residuum {

	namespace mesh {

		struct Data {

			// dimension
			Index parametricDim;
			Index spatialDim;

			// element/basis type
			ElementFamily elementFamily = ElementFamily::Quad;
			BasisType basisType = BasisType::Lagrange;

			// order
			std::vector<Index> basisOrder;

			// count
			Index numNodes;
			Index numElements;
			Index facesPerElement;
			Index nodesPerElement;

			// mesh data
			std::vector<Real> xyz;
			std::vector<Index> ien;
			std::vector<Int> rng;
			std::vector<Real> C;
		
		}; // struct Data
	
	} // namespace mesh

} // namespace residuum

#endif
