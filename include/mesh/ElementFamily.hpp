#ifndef RESIDUUM_MESH_ELEMENTFAMILY_HPP
#define RESIDUUM_MESH_ELEMENTFAMILY_HPP

namespace residuum {
	namespace mesh {

		enum class ElementFamily {
			Quad,
			Tri,
			Hex,
			Tet,
			Wedge
		}; // enum class ElementFamily

	} // namespace mesh
} // namespace residuum

#endif
