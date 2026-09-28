#ifndef RESIDUUM_FEM_BOUNDARY_BOUNDARYCATEGORY_HPP
#define RESIDUUM_FEM_BOUNDARY_BOUNDARYCATEGORY_HPP

namespace residuum {
	namespace fem {
		namespace boundary {

			enum class BCCategory {
				None,
				Essential,
				Natural
			}; // enum class BCCategory

		} // namespace boundary
	} // namespace fem
} // namespace residuum

#endif
