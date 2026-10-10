#ifndef RESIDUUM_APPLICATION_HEATEQ_CONFIG_MONITORCONFIG_HPP
#define RESIDUUM_APPLICATION_HEATEQ_CONFIG_MONITORCONFIG_HPP

#include <string>
#include <vector>

#include "fem/quantity/Reduction.hpp"

#include "solver/config/EvaluationMethod.hpp"
#include "solver/config/MonitorOutputConfig.hpp"
#include "solver/config/MonitorTermConfig.hpp"

namespace residuum {
	namespace application {
		namespace heateq {
			namespace config {

				struct MonitorConfig {

					enum class Quantity {
						HeatFlux
					}; // enum class Quantity

					std::string name;
					Quantity quantity;

					// only meaningful when quantity == HeatFlux
					solver::config::EvaluationMethod fluxMethod = solver::config::EvaluationMethod::Weak;

					fem::quantity::Reduction reduction;

					std::vector<solver::config::MonitorTermConfig> terms;

					solver::config::MonitorOutputConfig output;

				}; // struct MonitorConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace residuum

#endif
