#ifndef RESIDUUM_MESH_PARTITION_MESHPARTITIONER_HPP
#define RESIDUUM_MESH_PARTITION_MESHPARTITIONER_HPP

#include "core/Types.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/partition/PartitionData.hpp"
#include "mesh/partition/PartitionMethod.hpp"

namespace residuum {
	namespace mesh {
		namespace partition {

			class MeshPartitioner {
			public:
				PartitionData partition(const mesh::Mesh& mesh, Index numParts, PartitionMethod method = PartitionMethod::RCB);
			}; // class MeshPartitoner

		} // namespace partition
	} // namespace mesh
} // namespace residuum

#endif
