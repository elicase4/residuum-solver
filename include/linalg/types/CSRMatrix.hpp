#ifndef RESIDUUM_LINALG_TYPES_CSRMATRIX_HPP
#define RESIDUUM_LINALG_TYPES_CSRMATRIX_HPP

#include <memory>
#include <stdexcept>
#include <string>

#include "core/Types.hpp"

namespace residuum {
	namespace linalg {
		namespace types {
			
			template<typename T, typename BackendT>
			class CSRMatrix {
			public:

				using value_type = T;
				using backend_type = BackendT;

				CSRMatrix(Index nRows, Index nCols) : nRows_(nRows), nCols_(nCols), rowPtr_(BackendT::template alloc<Index>(nRows + 1)) {};
				
				// Move-only
				CSRMatrix(const CSRMatrix&) = delete;
				CSRMatrix& operator=(const CSRMatrix&) = delete;
				CSRMatrix(CSRMatrix&&) noexcept = default;
				CSRMatrix& operator=(CSRMatrix&&) noexcept = default;

				/* move to backend specialization */
				// Size
				Index nRows() const { return nRows_; }
				Index nCols() const { return nCols_; }

				// Resize
				void resize(Index nnz){
					nNnz_ = nnz;
					colIdx_ = BackendT::template alloc<Index>(nnz);
					data_ = BackendT::template alloc<T>(nnz);
				}

				// Access
				Index* rowPtr() { return rowPtr_.get(); }
				const Index* rowPtr() const { return rowPtr_.get(); }
				Index* colIdx() { return colIdx_.get(); }
				const Index* colIdx() const { return colIdx_.get(); }
				T* data() { return data_.get(); }
				const T* data() const { return data_.get(); }
				
				Index getDataIndex(Index i, Index j) const {

					Index start = rowPtr_.get()[i];
					Index end = rowPtr_.get()[i+1];

					const Index* cols = colIdx_.get();

					// binary search
					while (start < end){
						
						Index mid = start + (end - start) / 2;

						if (cols[mid] == j){
							return mid;
						} else if (cols[mid] < j){
							start = mid + 1;
						} else {
							end = mid;
						}

					}

					throw std::out_of_range("Index does not correspond to non-zero entry.");
					
				}

				// Zero-out data
				void zero(){
					BackendT::template zero<T>(data_.get(), nNnz_);
				}

			private:
				Index nRows_, nCols_, nNnz_;
				
				/* move to backend specialization */
				typename BackendT::template Ptr<Index> rowPtr_; // nRows + 1
				typename BackendT::template Ptr<Index> colIdx_; // nnz
				typename BackendT::template Ptr<T> data_; // nnz

			}; // class CSRMatrix
			
		} // namespace types
	} // namespace linalg
} // namespace residuum

#include "linalg/types/backend/CPU.hpp"
//#include "linalg/types/backend/CUDA.hpp"
		
#endif
