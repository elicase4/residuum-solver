#ifndef RESIDUUM_APPLICATION_HEATEQ_CONFIG_OUTPUTFIELDCONFIG_HPP
#define RESIDUUM_APPLICATION_HEATEQ_CONFIG_OUTPUTFIELDCONFIG_HPP

#include <vector>

namespace residuum {
	namespace application {
		namespace heateq {
			namespace config {

				enum class OutputFieldQuantity {
					HeatFlux,
					Conductivity,
					SpecificHeat,
					Residual
				}; // enum class OutputFieldQuantity

				struct OutputFieldConfig {

					std::vector<OutputFieldQuantity> fields;

				}; // struct OutputFieldConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace residuum

#endif
