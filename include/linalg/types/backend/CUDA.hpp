#ifndef RESIDUUM_LINALG_TYPES_BACKEND_CUDA_HPP
#define RESIDUUM_LINALG_TYPES_BACKEND_CUDA_HPP

#include <cuda_runtime.h>
#include <memory>
#include <stdexcept>

#include "core/Types.hpp"
#include "linalg/types/backend/cuda/CopyKind.hpp"

namespace residuum {
	namespace linalg {
		namespace types {
			namespace backend {
				
				template<typename T>
				struct CudaDeleter {
					void operator()(T* ptr) const noexcept {
						if (ptr) cudaFree(ptr);
					}
				}; // struct CudaDeleter

				class CUDA {
				public:
					
					template<typename T>
					using Ptr = std::unique_ptr<T, CudaDeleter<T>>;

					template<typename T>
					static Ptr<T> alloc(Index n){
						T* ptr = nullptr;
						cudaError_t err = cudaMalloc(&ptr, n*sizeof(T));
						if (err != cudaSuccess){
							throw std::runtime_error("cudaMalloc failure");
						}
						return Ptr<T>(ptr);
					}

					template<typename T>
					static void copy(T* dst, const T* src, Index n, CopyKind kind){
						
						cudaMemcpyKind CUDAKind;
						
						switch(kind) {
							case (CopyKind::HostToDevice):
								CUDAKind = cudaMemcpyHostToDevice;
								break;
							case (CopyKind::DeviceToHost):
								CUDAKind = cudaMemcpyDeviceToHost;
								break;
							case (CopyKind::DeviceToDevice):
								CUDAKind = cudaMemcpyDeviceToDevice;
								break;
							default:
								CUDAKind = cudaMemcpyDefault;
						}

						cudaMemcpy(dst, src, n*sizeof(T), CUDAkind);
					}
					
					template<typename T>
					static void zero(T* ptr, Index n){
						cudaMemset(ptr, 0, n*sizeof(T));
					}
					
					/*
					template<typename T>
					static void set(T* ptr, Index n, T value){
						cudaMemset(ptr, value, n*sizeof(T));
					}
					*/

				}; // class CUDA

			} // namespace backend
		} // namespace types
	} // namespace linalg
} // namespace residuum

#endif
