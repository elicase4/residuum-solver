#ifndef RESIDUUM_APPLICATION_HEATEQ_CONFIG_CONDUCTIVITYCONFIG_HPP
#define RESIDUUM_APPLICATION_HEATEQ_CONFIG_CONDUCTIVITYCONFIG_HPP

#include <string>
#include <vector>

#include "core/Types.hpp"

namespace residuum {
	namespace application {
		namespace heateq {
			namespace config {

				struct ConductivityConfig {

					enum class Type {
						Default,
						Constant,
						Anisotropic,
						TemperatureDependentIsotropic,
						TemperatureDependentAnisotropic
					}; // enum class Type

					Type type;

					Real value; // Constant

					std::vector<std::vector<Real>> tensor; // Anisotropic, TemperatureDependentAnisotropic

					std::string valueExpression;    // TemperatureDependentIsotropic/Anisotropic: k(T)
					std::string gradientExpression; // TemperatureDependentIsotropic/Anisotropic: dk/dT(T)

					std::string unit; // required, SI unit label

				}; // struct ConductivityConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace residuum

#endif
