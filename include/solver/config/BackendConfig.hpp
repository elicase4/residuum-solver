#ifndef RESIDUUM_SOLVER_CONFIG_BACKENDCONFIG_HPP
#define RESIDUUM_SOLVER_CONFIG_BACKENDCONFIG_HPP

namespace residuum {
	namespace solver {
		namespace config {

			struct BackendConfig {

				enum class Type {
					CPU,
					CUDA
				}; // enum class Type

				Type type = Type::CPU;

			}; // struct BackendConfig

		} // namespace config
	} // namespace solver
} // namespace residuum

#endif
