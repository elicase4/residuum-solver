#ifndef RESIDUUM_IO_FIELDIO_NODALVALUESOURCEADAPTER_HPP
#define RESIDUUM_IO_FIELDIO_NODALVALUESOURCEADAPTER_HPP

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include "core/Types.hpp"
#include "fem/evaluator/EvalNodalData.hpp"
#include "mesh/Mesh.hpp"

namespace residuum {
	namespace io {
		namespace fieldio {

			template<fem::evaluator::EvalNodalData SourceT>
			class NodalValueSourceAdapter {
			public:

				static constexpr Index NumComponents = SourceT::NumComponents;

				template<typename... Args>
				NodalValueSourceAdapter(const mesh::Mesh& mesh, Args&&... args) : source_(mesh, std::forward<Args>(args)...), spatialDim_(mesh.data.spatialDim), sortedNodeIndices_(mesh.data.numNodes) {

					// copy xyz into flatCoords
					flatCoords_ = mesh.data.xyz;

					// fill sortedNodeIndices
					std::iota(sortedNodeIndices_.begin(), sortedNodeIndices_.end(), 0);

					// is index a less than index b?
					auto compareNodes = [&](Index a, Index b) {
						
						for (Index j = 0; j < spatialDim_; ++j) {
							if (flatCoords_[a*spatialDim_ + j] > flatCoords_[b*spatialDim_ + j]) {
								return false;
							} else if (flatCoords_[b*spatialDim_ + j] > flatCoords_[a*spatialDim_ + j]) {
								return true;
							}
						}

						return false;
					
					}; // compareNodes

					std::sort(sortedNodeIndices_.begin(), sortedNodeIndices_.end(), compareNodes);

				}

				void operator()(Real, const Real* x, Real* out) const {
				
					auto compareNodeWithLessIndex = [&](Index mid, const Real* x) {
						
						for (Index j = 0; j < spatialDim_; ++j) {
							if (flatCoords_[mid*spatialDim_ + j] < x[j]){
								return true;
							} else if (flatCoords_[mid*spatialDim_ + j] > x[j]) {
								return false;
							}
						}
						
						return false;
					
					}; // compareNodeWithIndex

					// does the nodeID a match the inlet coordinate?
					auto isCoordinateAtNode = [&](Index nodeID) {
						
						for (Index j = 0; j < spatialDim_; ++j) {
							if (flatCoords_[nodeID*spatialDim_ + j] != x[j]){
								return false;
							}
						}
						
						return true;
					
					}; // isCoordinateAtNode

					// binary search
					auto it = std::lower_bound(sortedNodeIndices_.begin(), sortedNodeIndices_.end(), x, compareNodeWithLessIndex);

					if (it != sortedNodeIndices_.end() && isCoordinateAtNode(*it)){
						source_.eval(*it, out);
					} else {
						throw std::runtime_error("Nodal Source: Node coordinate could not be located.");
					}

				}

			private:

				SourceT source_;
				Index spatialDim_;
				std::vector<Real> flatCoords_;       // size numNodes*spatialDim_, node i's coords at [i*spatialDim_, (i+1)*spatialDim_)
				std::vector<Index> sortedNodeIndices_; // node indices into flatCoords_, sorted lexicographically by coordinate

			}; // class NodalValueSourceAdapter

		} // namespace fieldio
	} // namespace io
} // namespace residuum

#endif
