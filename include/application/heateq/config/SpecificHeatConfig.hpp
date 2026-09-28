#ifndef RESIDUUM_APPLICATION_HEATEQ_CONFIG_SPECIFICHEATCONFIG_HPP
#define RESIDUUM_APPLICATION_HEATEQ_CONFIG_SPECIFICHEATCONFIG_HPP

#include <string>

#include "core/Types.hpp"

namespace residuum {
	namespace application {
		namespace heateq {
			namespace config {

				struct SpecificHeatConfig {

					enum class Type {
						Default,
						Constant,
						TemperatureDependent
					}; // enum class Type

					Type type;

					Real value; // Constant

					std::string valueExpression;    // TemperatureDependent: cp(T)
					std::string gradientExpression; // TemperatureDependent: dcp/dT(T)

					std::string unit; // required, SI unit label

				}; // struct SpecificHeatConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace residuum

#endif
