#ifndef RESIDUUM_LINALG_TYPES_BACKEND_CUDA_COPYKIND_HPP
#define RESIDUUM_LINALG_TYPES_BACKEND_CUDA_COPYKIND_HPP

namespace residuum {
	namespace linalg {
		namespace types {
			namespace backend {
				
				enum class CopyKind {
					HostToHost,
					HostToDevice,
					DeviceToHost,
					DeviceToDevice,
				}; // class CopyKind
				
			} // namespace backend
		} // namespace types
	} // namespace linalg
} // namespace residuum

#endif
