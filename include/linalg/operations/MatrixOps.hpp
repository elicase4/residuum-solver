#ifndef RESIDUUM_LINALG_OPERATIONS_MATRIXOPS_HPP
#define RESIDUUM_LINALG_OPERATIONS_MATRIXOPS_HPP

#include "linalg/types/Vector.hpp"

namespace residuum {
	namespace linalg {
		namespace operations {

			template<typename MatrixT, typename VectorT>
			void matvec(const MatrixT& A, const VectorT& x, VectorT& y);

			template<typename MatrixT, typename VectorT>
			void lump(const MatrixT& A, VectorT& diag);

		} // namespace operations
	} // namespace linalg
} // namespace residuum


#include "linalg/operations/backend/cpu/MatrixOps.tpp"
//#include "linalg/operations/backend/cuda/MatrixOps.tpp"

#endif
