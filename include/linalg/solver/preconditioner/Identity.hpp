#ifndef RESIDUUM_LINALG_SOLVER_PRECONDITIONER_IDENTITY_HPP
#define RESIDUUM_LINALG_SOLVER_PRECONDITIONER_IDENTITY_HPP

#include "linalg/operations/VectorOps.hpp"

namespace residuum {
	namespace linalg {
		namespace solver {
			namespace preconditioner {

				template<typename VectorT>
				class Identity {
				public:
					
					void apply(const VectorT& r, VectorT& z){
						operations::copy(r, z); // z = r
					}

					// no state to refresh against the current operator
					template<typename OperatorT>
					void update(const OperatorT&) {}

					Index flopsPerApply() const {
						return 0;
					}

				}; // class Identity

			} // namespaace preconditioner
		} // namespaace solver
	} // namespaace linalg
} // namespaace residuum

#endif
