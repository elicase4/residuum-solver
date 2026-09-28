#ifndef RESIDUUM_LINALG_TYPES_BACKEND_CPU_HPP
#define RESIDUUM_LINALG_TYPES_BACKEND_CPU_HPP

#include <memory>
#include <algorithm>

#include "core/Types.hpp"
#include "linalg/types/backend/cuda/CopyKind.hpp"

namespace residuum {
	namespace linalg {
		namespace types {
			namespace backend {
				
				class CPU {
				public:
					template<typename T>
					using Ptr = std::unique_ptr<T[]>;

					template<typename T>
					static Ptr<T> alloc(Index n){
						return std::make_unique<T[]>(n);
					}

					template<typename T>
					static void copy(T* dst, const T* src, Index n, CopyKind){
						std::copy(src, src + n, dst);
					}

					template<typename T>
					static void zero(T* ptr, Index n){
						std::fill(ptr, ptr + n, T(0));
					}

					template<typename T>
					static void set(T* ptr, Index n, T value){
						std::fill(ptr, ptr + n, T(value));
					}

				}; // class CPU

			} // namespace backend
		} // namespace types
	} // namespace linalg
} // namespace residuum

#endif
