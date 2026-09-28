#ifndef RESIDUUM_LINALG_TYPES_VECTOR_HPP
#define RESIDUUM_LINALG_TYPES_VECTOR_HPP

#include "core/Types.hpp"

namespace residuum {
	namespace linalg {
		namespace types {

			template<typename T, typename BackendT>
			class Vector {
			public:
				
				using value_type = T;
				using backend_type = BackendT;

				explicit Vector(Index size) : size_(size), data_(BackendT::template alloc<T>(size)) {};
				
				// Move only operations
				Vector(const Vector&) = delete;
				Vector& operator=(const Vector&) = delete;
				Vector(Vector&&) = default;
				Vector& operator=(Vector&&) = default;

				// Size
				Index size() const { return size_; }

				// Access
				T* data() { return data_.get(); }
				const T* data() const { return data_.get(); }

				// Zero-out data
				void zero(){
					BackendT::template zero<T>(data_.get(), size_);
				}

				void set(T value){
					BackendT::template set<T>(data_.get(), size_, value);
				}
			
			private:
				Index size_;
				typename BackendT::template Ptr<T> data_;

			}; // class Vector
		
		} // namespace types
	} // namespace linalg
} // namespace residuum

#include "linalg/types/backend/CPU.hpp"
//#include "linalg/types/backend/CUDA.hpp"

#endif
