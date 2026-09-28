#ifndef RESIDUUM_MESH_MESH_HPP
#define RESIDUUM_MESH_MESH_HPP

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "core/Types.hpp"
#include "mesh/Data.hpp"

namespace residuum {
	
	namespace mesh {
		
		class Mesh {
		public:
				
			Mesh() = default;

			Real* getNodeCoord(Index nodeID) { return &data.xyz[nodeID * data.spatialDim]; }
			const Real* getNodeCoord(Index nodeID) const { return &data.xyz[nodeID * data.spatialDim]; }
			
			Index* getElementNodes(Index elemID) { return &data.ien[elemID * data.nodesPerElement]; }
			const Index* getElementNodes(Index elemID) const { return &data.ien[elemID * data.nodesPerElement]; }

			Int* getBoundaryTag(Index elemID) { return &data.rng[elemID * data.facesPerElement]; }
			const Int* getBoundaryTag(Index elemID) const { return &data.rng[elemID * data.facesPerElement]; }
			bool isOnBoundary(Index elemID, Index faceID) const { return data.rng[elemID * data.facesPerElement + faceID] >= 0; }
			bool isOnBoundaryTag(Index elemID, Index faceID, Int tag){ return data.rng[elemID * data.facesPerElement + faceID] == tag; }
			
			//Index extractionMatrixSize() const{};
			//Real* getExtractionOperator(Index elemID);
			//const Real* getExtractionOperator(Index elemID) const;
			
			void clear() {
				data.parametricDim = 0; data.spatialDim = 0; data.numNodes = 0; data.numElements = 0; data.nodesPerElement = 0; data.facesPerElement = 0;
				data.xyz.clear(); data.ien.clear(); data.rng.clear(); data.C.clear(); data.basisOrder.clear();
			}

			bool isValid() const {
				
				if (data.xyz.size() != data.numNodes * data.spatialDim) {
					return false;
				}
				if (data.ien.size() != data.numElements * data.nodesPerElement) {
					return false;
				}
				if (data.rng.size() != data.numElements * data.facesPerElement) {
					return false;
				}
				if (data.parametricDim == 0 || data.spatialDim == 0) {
					return false;
				}
				for (Index p : data.basisOrder){
					if (p == 0){
						return false;
					}
				}
				if (data.nodesPerElement != computeNodesPerElement(data.basisOrder)) {
					return false;
				}
				
				return true;
			}

			void print(std::ostream& os = std::cout) const;
			
			bool isIGA() const { return !data.C.empty(); }
		
			Data data;
		
		protected:

			static Index computeNodesPerElement(std::vector<Index> basisOrder) {
				Index npe = 1;
				for (Index p : basisOrder){
					npe *= (p + 1);
				}
				return npe;
			}
			
		}; // class Mesh
	
	} // namespace mesh

} // namespace residuum

#endif
