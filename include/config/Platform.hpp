#ifndef RESIDUUM_CONFIG_PLATFORM_HPP
#define RESIDUUM_CONFIG_PLATFORM_HPP

// device accessible functions
#ifdef __CUDACC__
	#define PDE_HOST __host__
	#define PDE_DEVICE __device__
	#define PDE_INLINE __forceinline__
#else
	#define PDE_HOST
	#define PDE_DEVICE
	#define PDE_INLINE inline
#endif

// hints to the optimizer that a branch (e.g. an order/count dispatch already bounds-checked
// at a system boundary such as fem::dispatch::validateDiscretizationLimits) cannot be reached
#define PDE_UNREACHABLE() __builtin_unreachable()

#endif
