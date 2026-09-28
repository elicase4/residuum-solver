#ifndef RESIDUUM_LINALG_TYPES_MATRIX_HPP
#define RESIDUUM_LINALG_TYPES_MATRIX_HPP

#include "core/Types.hpp"

namespace residuum {
	namespace linalg {
		namespace types {

			template<typename T, typename BackendT>
			class Matrix {
			public:

				using value_type = T;
				using backend_type = BackendT;

				explicit Matrix(Index nRows, Index nCols) : nRows_(nRows), nCols_(nCols), data_(BackendT::template alloc<T>(nRows * nCols)) {};
				
				// Move only operations
				Matrix(const Matrix&) = delete;
				Matrix& operator=(const Matrix&) = delete;
				Matrix(Matrix&&) = default;
				Matrix& operator=(Matrix&&) = default;

				// Size
				Index nRows() const { return nRows_; }
				Index nCols() const { return nCols_; }

				// Access
				T* data() { return data_.get(); }
				const T* data() const { return data_.get(); }
				
				// Zero-out data
				void zero(){
					BackendT::template zero<T>(data_.get(), nRows_ * nCols_);
				}

			private:
				Index nRows_, nCols_;
				typename BackendT::template Ptr<T> data_;

			}; // class Matrix
		
			} // namespace types
		} // namespace linalg
	} // namespace residuum

#include "linalg/types/backend/CPU.hpp"
//#include "linalg/types/backend/CUDA.hpp"

#endif
