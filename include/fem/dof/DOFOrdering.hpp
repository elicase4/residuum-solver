#ifndef RESIDUUM_FEM_DOF_DOFORDERING_HPP
#define RESIDUUM_FEM_DOF_DOFORDERING_HPP

namespace residuum {
	namespace fem {
		namespace dof {

			enum class DOFOrdering {
				Interleaved, // node-major
				Block // field-major
			}; // enum class DOFOrdering

		} // namespace dof
	} // namespace fem
} // namespace residuum

#endif
