#include <cassert>

namespace residuum::linalg::operations {

	template<typename MatrixT, typename VectorT>
	void matvec(const MatrixT& A, const VectorT& x, VectorT& y){

		assert(A.nCols() == x.size());
		assert(A.nRows() == y.size());

		for (Index i = 0; i < A.nRows(); ++i){

			typename VectorT::value_type sum = 0;

			for (Index p = A.rowPtr()[i]; p < A.rowPtr()[i+1]; ++p){
				sum += A.data()[p] * x.data()[A.colIdx()[p]];
			}

			y.data()[i] = sum;

		}

	}

	template<typename MatrixT, typename VectorT>
	void lump(const MatrixT& A, VectorT& diag){

		assert(A.nCols() == A.nCols());
		assert(A.nRows() == diag.size());

		for (Index i = 0; i < A.nRows(); ++i){

			typename VectorT::value_type sum = 0;

			for (Index p = A.rowPtr()[i]; p < A.rowPtr()[i+1]; ++p){
				sum += A.data()[p];
			}

			diag.data()[i] = sum;
		
		}

	}

} // namespace residuum::linalg::operations
