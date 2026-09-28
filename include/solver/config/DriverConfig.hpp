#ifndef RESIDUUM_SOLVER_CONFIG_DRIVERCONFIG_HPP
#define RESIDUUM_SOLVER_CONFIG_DRIVERCONFIG_HPP

namespace residuum {
	namespace solver {
		namespace config {

			struct DriverConfig {

				enum class Type {
					Steady,
					Transient,
					PseudoTransient // steady, driven by pseudo-time-marching, not yet implemented
				}; // enum class Type

				Type type = Type::Steady;

			}; // struct DriverConfig

		} // namespace config
	} // namespace solver
} // namespace residuum

#endif
