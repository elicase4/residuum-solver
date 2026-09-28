#ifndef RESIDUUM_IO_FIELDIO_NODALVALUESOURCEADAPTER_HPP
#define RESIDUUM_IO_FIELDIO_NODALVALUESOURCEADAPTER_HPP

#include <map>
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
				NodalValueSourceAdapter(const mesh::Mesh& mesh, Args&&... args) : source_(mesh, std::forward<Args>(args)...), spatialDim_(mesh.data.spatialDim) {

					for (Index nodeID = 0; nodeID < mesh.data.numNodes; ++nodeID) {
						const Real* xyz = mesh.getNodeCoord(nodeID);
						coordToNode_[std::vector<Real>(xyz, xyz + spatialDim_)] = nodeID;
					}

				}

				void operator()(Real time, const Real* x, Real* out) const {

					(void) time;

					auto it = coordToNode_.find(std::vector<Real>(x, x + spatialDim_));

					if (it == coordToNode_.end()) {
						throw std::runtime_error("NodalValueSourceAdapter::operator(): coordinate not found in mesh, x must come from mesh.getNodeCoord() on the mesh this adapter was constructed against");
					}

					source_.eval(it->second, out);

				}

			private:

				SourceT source_;
				Index spatialDim_;
				std::map<std::vector<Real>, Index> coordToNode_;

			}; // class NodalValueSourceAdapter

		} // namespace fieldio
	} // namespace io
} // namespace residuum

#endif
